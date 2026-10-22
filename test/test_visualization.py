"""Check scientific field layout without allocating a full-volume browser dataset."""
import tempfile
import unittest
from pathlib import Path

from tools.export_demo_data import grain_palette, read_vtk


class ExportGeometry(unittest.TestCase):
    def test_streamed_planes_and_exterior_faces(self):
        nx, ny, nz = 4, 5, 3
        field = lambda i, j, k: i + 10*j + 100*k
        ids = [field(i, j, k) for k in range(nz) for j in range(ny) for i in range(nx)]
        with tempfile.TemporaryDirectory() as directory:
            vtk = Path(directory) / "field.vtk"
            body = "DIMENSIONS 5 6 4\nSPACING 0.0005 0.0005 0.0005\n"
            for name, values in (("grain_id", ids), ("cell_state", [1]*len(ids)),
                                 ("temperature_K", [1500 + i/10 for i in ids])):
                body += f"SCALARS {name} double 1\nLOOKUP_TABLE default\n"
                body += "\n".join(map(str, values)) + "\n"
            vtk.write_text(body)
            dimensions, spacing, planes, result = read_vtk(vtk, 2)
        self.assertEqual(dimensions, [nx, ny, nz])
        self.assertEqual(spacing, .0005)
        self.assertEqual(planes, [0, 2, 4])
        selected = [field(i,j,k) for k in range(nz) for j in planes for i in range(nx)]
        self.assertEqual(result["grain_id"], selected)
        self.assertEqual(result["temperature_K"], [1500 + i/10 for i in selected])
        self.assertEqual(result["surfaces"]["top"], [field(i,j,nz-1) for j in range(ny) for i in range(nx)])
        self.assertEqual(result["surfaces"]["side"], [field(nx-1,j,k) for k in range(nz) for j in range(ny)])
        self.assertEqual(result["assigned_fraction"], (len(ids)-1)/len(ids))

    def test_palette_keeps_unassigned_distinct(self):
        palette = grain_palette(100)
        self.assertEqual(len(palette), 101)
        self.assertEqual(len(set(palette)), 101)
        self.assertNotIn(palette[0], palette[1:])


if __name__ == "__main__":
    unittest.main()
