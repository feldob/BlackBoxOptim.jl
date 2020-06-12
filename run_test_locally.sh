docker run -it --rm -e JULIA_NUM_THREADS=4 -v "$PWD":/usr/myapp -w /usr/myapp testbbo julia --color=yes test/runtests.jl
