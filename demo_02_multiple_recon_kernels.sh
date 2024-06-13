# run with: `bash demo_02_multiple_recon_kernels.sh`
#
# Shell script example of how to reproduce sets of simulations with changing parameters using a config file, here 'configs/multiple_recon_kernels.toml'
# See Usage documentation for more: <https://pediatriciqphantoms.readthedocs.io/en/latest/usage.html#command-line-interface>
config='configs/multiple_recon_kernels.toml'

make_phantoms $config