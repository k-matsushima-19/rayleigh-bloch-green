# rayleigh-bloch-green
C++ code for computing Rayleigh-Bloch waves swapping between physical and unphysical Riemann sheets

# How to run on Docker
```bash
git clone https://github.com/k-matsushima-19/rayleigh-bloch-green.git
cd rayleigh-bloch-green
mkdir output
docker build . -t rayleigh-bloch-green
docker run -it --rm --mount type=bind,src=$PWD/output,dst=/app/output rayleigh-bloch-green
```