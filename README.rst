High-Throughput COVID Fingerprinting Pipeline
=============================================

Starting with fp2.py code from @ianfoster.

Setting up your Theta env.
--------------------------

Step 1: Setup conda env

>>> module load miniconda-3/latest
>>> conda create -p /projects//candle_aesp/yadu/candle_inferpy3.7 --clone $CONDA_PREFIX
>>> conda activate /projects/candle_aesp/yadu/candle_inferpy3.7

Step 2: Install all the required packaged. Most packages already come baked in with the base conda

>>> git clone https://github.com/Parsl/parsl.git
>>> cd parsl
>>> pip install .

>>> pip install conda-pack

Running the pipeline
--------------------

