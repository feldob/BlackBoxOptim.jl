FROM julia:latest

RUN julia -e 'using Pkg; Pkg.add("Distributions")'
RUN julia -e 'using Pkg; Pkg.add("StatsBase")'
RUN julia -e 'using Pkg; Pkg.add("Compat")'
RUN julia -e 'using Pkg; Pkg.add("SpatialIndexing")'
RUN julia -e 'using Pkg; Pkg.add("CPUTime")'
RUN julia -e 'using Pkg; Pkg.add("HTTP")'
RUN julia -e 'using Pkg; Pkg.add("JSON")'

