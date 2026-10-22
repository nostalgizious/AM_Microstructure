searchdir(path,key) = filter(x->occursin(key,x), readdir(path)) #Search for a file in a directory
wait_for_key(prompt) = (print(stdout, prompt); read(stdin, 1); nothing)
# wait_for_key(prompt) = print(prompt)