#lib/python3.6/site-packages/graph_tool/inference/mcmc.py

python_v=`python -V | cut -f 1,2 -d. | awk '{print $2}'`
cwd=`cwd`

pushd $CONDA_PREFIX/lib/python${python_v}/site-packages/graph_tool/inference
patch -p0 < ${cwd}/patch.mcmc.txt

popd

python conftest.py
