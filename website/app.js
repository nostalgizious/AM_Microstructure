"use strict";
(() => {
  const data = window.AM_DEMO;
  if (!data) return;
  const [nx, ny, nz] = data.dimensions;
  const sliceIndices=data.slice_indices, ns=sliceIndices.length;
  const palette = data.palette;
  const color = id => palette[id];
  const at = (x,y,z) => x + nx*(y + ns*z);
  const volume = document.getElementById("volume");
  const vc = volume.getContext("2d");
  const final = data.frames.at(-1);
  const rawProject = (x,y,z) => [x + y*.8, x*.15 - y*.5 - z];
  const corners = [0,nx].flatMap(x => [0,ny].flatMap(y => [0,nz].map(z => rawProject(x,y,z))));
  const minX = Math.min(...corners.map(p=>p[0])), maxX = Math.max(...corners.map(p=>p[0]));
  const minY = Math.min(...corners.map(p=>p[1])), maxY = Math.max(...corners.map(p=>p[1]));
  const scale = Math.min((volume.width-150)/(maxX-minX),(volume.height-180)/(maxY-minY));
  const offsetX=(volume.width-scale*(maxX-minX))/2, offsetY=(volume.height-120-scale*(maxY-minY))/2;
  const project = (x,y,z) => {const p=rawProject(x,y,z);return [offsetX+(p[0]-minX)*scale,offsetY+(p[1]-minY)*scale]};
  const quad = (points, fill, shade=0) => {
    vc.beginPath(); points.forEach((p,i) => i ? vc.lineTo(...p) : vc.moveTo(...p)); vc.closePath(); vc.fillStyle=fill; vc.fill();
    if(shade){vc.fillStyle=`rgba(16,40,24,${shade})`;vc.fill()}
    // A subpixel stroke closes canvas seams without adding a voxel grid.
    vc.strokeStyle=fill;vc.lineWidth=.3;vc.stroke();
  };
  vc.clearRect(0,0,volume.width,volume.height);
  // Only the visible boundary faces of the real 3D voxel grid are rendered.
  for(let y=ny-1;y>=0;y--) for(let x=0;x<nx;x++) quad([project(x,y,nz),project(x+1,y,nz),project(x+1,y+1,nz),project(x,y+1,nz)],color(final.surfaces.top[x+nx*y]));
  for(let y=ny-1;y>=0;y--) for(let z=0;z<nz;z++) quad([project(nx,y,z),project(nx,y+1,z),project(nx,y+1,z+1),project(nx,y,z+1)],color(final.surfaces.side[y+ny*z]),.2);
  for(let z=0;z<nz;z++) for(let x=0;x<nx;x++) quad([project(x,0,z),project(x+1,0,z),project(x+1,0,z+1),project(x,0,z+1)],color(final.grain_id[at(x,0,z)]),.04);
  const bar=20/(data.spacing_mm*1000)*scale;
  vc.strokeStyle="#4b6744";vc.lineWidth=2;vc.beginPath();vc.moveTo(90,810);vc.lineTo(90+bar,810);vc.stroke();vc.fillStyle="#4b6744";vc.font="22px monospace";vc.fillText("20 µm",90,848);vc.fillText("z ↑",25,420);vc.fillText("x →",850,848);
  document.documentElement.classList.add("js-ready");
  const slice = document.getElementById("slice");
  // One canvas pixel per actual cell avoids antialiased seams between cells.
  slice.width=nx;slice.height=nz;slice.parentElement.style.aspectRatio=`${nx}/${nz}`;
  const sc = slice.getContext("2d");
  const yInput=document.getElementById("slice-position"), field=document.getElementById("field");
  yInput.max=String(ns-1);yInput.value=String(sliceIndices.indexOf(data.summary.slice_y_index-1));
  let stage=data.frames.length-1;
  const tempColor=t=>{const v=Math.max(0,Math.min(1,(t-1500)/300));return `hsl(${(1-v)*220},65%,${35+v*25}%)`};
  function draw(){
    const frame=data.frames[stage], y=Number(yInput.value), cw=slice.width/nx,ch=slice.height/nz;
    for(let z=0;z<nz;z++)for(let x=0;x<nx;x++){
      const index=at(x,y,z);sc.fillStyle=field.value==="grain"?color(frame.grain_id[index]):tempColor(frame.temperature_K[index]);sc.fillRect(x*cw,(nz-z-1)*ch,cw,ch);
    }
    const fraction=frame.assigned_fraction;
    document.getElementById("viewer-time").textContent=frame.time_ms.toFixed(2)+" ms";
    document.getElementById("assigned").textContent=(100*fraction).toFixed(1)+"% assigned";
    document.getElementById("slice-label").textContent=((sliceIndices[y]+.5)*data.spacing_mm*1000).toFixed(2)+" µm";
    document.getElementById("field-legend").textContent=field.value==="grain"?"Categorical grain colors · dark = no assigned grain":"Temperature: blue 1500 K → orange/red 1800 K";
    slice.setAttribute("aria-label",`${frame.name} stage at ${frame.time_ms} milliseconds, physical y slice ${sliceIndices[y]+1} of ${ny}, ${field.value} field. ${(fraction*100).toFixed(1)} percent of cells have grain IDs.`);
  }
  document.querySelectorAll("[data-stage]").forEach(button=>button.addEventListener("click",()=>{stage=Number(button.dataset.stage);document.querySelectorAll("[data-stage]").forEach(b=>b.setAttribute("aria-pressed",String(b===button)));draw()}));
  yInput.addEventListener("input",draw);field.addEventListener("change",draw);draw();
  slice.addEventListener("pointermove",event=>{const r=slice.getBoundingClientRect();const x=Math.min(nx-1,Math.max(0,Math.floor((event.clientX-r.left)/r.width*nx))),z=Math.min(nz-1,Math.max(0,nz-1-Math.floor((event.clientY-r.top)/r.height*nz))),y=Number(yInput.value),i=at(x,y,z),frame=data.frames[stage];document.getElementById("cell-tooltip").textContent=`Cell (${x+1}, ${sliceIndices[y]+1}, ${z+1}) · Grain ${frame.grain_id[i] || "unassigned"} · ${frame.temperature_K[i].toFixed(1)} K`});
  const dialog=document.getElementById("figure-dialog");let opener;
  document.querySelectorAll("[data-zoom]").forEach(a=>a.addEventListener("click",e=>{if(!dialog.showModal)return;e.preventDefault();opener=a;document.getElementById("dialog-image").src=a.href;document.getElementById("dialog-image").alt=a.querySelector("img").alt;document.getElementById("dialog-caption").textContent=a.closest("figure").querySelector("figcaption").textContent;dialog.showModal()}));
  dialog.querySelector("button").addEventListener("click",()=>dialog.close());dialog.addEventListener("click",e=>{if(e.target===dialog)dialog.close()});dialog.addEventListener("close",()=>opener?.focus());
  document.getElementById("copy-command").addEventListener("click",async e=>{const text=document.getElementById("replay-command").textContent;try{await navigator.clipboard.writeText(text);e.target.textContent="Copied";setTimeout(()=>e.target.textContent="Copy command",1800)}catch{const range=document.createRange();range.selectNodeContents(document.getElementById("replay-command"));const selection=window.getSelection();selection.removeAllRanges();selection.addRange(range);e.target.textContent="Select & copy"}});
})();

// Playback uses an image exported from the real CA at every recorded timestep.
(() => {
  const data = window.AM_EVOLUTION;
  const slider = document.getElementById("evolution-slider");
  if (!data || !slider) return;
  const img = document.getElementById("evolution-frame");
  const play = document.getElementById("evolution-play");
  const speed = document.getElementById("evolution-speed");
  const last = data.frames.length - 1;
  const nearestFrame=time=>data.frames.reduce((best,frame,i)=>Math.abs(frame.time_ms-time)<Math.abs(data.frames[best].time_ms-time)?i:best,0);
  let index = nearestFrame(2.0), timer = null;
  const bufferedImages=new Map();
  const framePath=i=>`assets/demo/evolution/frame-${String(data.frames[i].step).padStart(3, "0")}.png`;
  function bufferFrames() {
    // A short look-ahead keeps playback smooth without decoding the full movie.
    for(let i=index;i<=Math.min(last,index+6);i++) {
      if(!bufferedImages.has(i)){const image=new Image();image.src=framePath(i);bufferedImages.set(i,image)}
    }
    for(const i of bufferedImages.keys())if(i<index-1 || i>index+6)bufferedImages.delete(i);
  }
  slider.max = String(last);
  document.querySelectorAll('[data-frame-time]').forEach(button=>button.dataset.index=String(nearestFrame(Number(button.dataset.frameTime))));

  function render() {
    const frame = data.frames[index];
    const path = framePath(index);
    img.src = path;
    img.alt = `Actual grain field at ${frame.time_ms.toFixed(3)} milliseconds, x-z slice at y = ${data.y_um} micrometres.`;
    slider.value = String(index);
    document.getElementById("evolution-time").textContent = `t = ${frame.time_ms.toFixed(3)} ms · step ${frame.step} / ${data.frames[last].step}`;
    const cycle=Math.max(0,data.heating_intervals_ms.findLastIndex(interval=>frame.time_ms>interval[0]+1e-8));
    const heating=frame.time_ms<=data.heating_intervals_ms[cycle][1]+1e-8;
    document.getElementById("evolution-phase").textContent = index===0 ? "Initial solid substrate" : `Pass ${cycle+1} · ${heating ? "heating and remelting" : "cooling and grain regrowth"}`;
    document.getElementById("evolution-statistics").textContent = `${frame.grains} grains intersect this section · ${(frame.assigned_fraction*100).toFixed(1)}% assigned`;
    document.getElementById("evolution-download").href = path;
    document.getElementById("evolution-prev").disabled = index === 0;
    document.getElementById("evolution-next").disabled = index === last;
    bufferFrames();
  }
  function pause() {
    clearTimeout(timer);
    timer = null;
    play.textContent = "Play";
    play.setAttribute("aria-pressed", "false");
  }
  function tick() {
    if (index === last) { pause(); return; }
    const next=bufferedImages.get(index+1);
    if(next && !next.complete){timer=setTimeout(tick,20);return}
    if(next && !next.naturalWidth){pause();return}
    index += 1;
    render();
    timer = setTimeout(tick, 1000 / Number(speed.value));
  }
  play.addEventListener("click", () => {
    if (timer !== null) { pause(); return; }
    if (index === last) { index = 0; render(); }
    play.textContent = "Pause";
    play.setAttribute("aria-pressed", "true");
    timer = setTimeout(tick, 1000 / Number(speed.value));
  });
  slider.addEventListener("input", () => { pause(); index = Number(slider.value); render(); });
  document.getElementById("evolution-prev").addEventListener("click", () => { pause(); index = Math.max(0, index - 1); render(); });
  document.getElementById("evolution-next").addEventListener("click", () => { pause(); index = Math.min(last, index + 1); render(); });
  document.querySelectorAll("[data-frame-time]").forEach(button => button.addEventListener("click", () => {
    pause(); index = Math.min(last, Number(button.dataset.index)); render();
  }));
  document.addEventListener("visibilitychange", () => { if (document.hidden) pause(); });
  render();
})();
