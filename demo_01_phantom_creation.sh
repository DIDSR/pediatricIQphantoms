# run with: `bash demo_01_phantom_creation.sh`
#
# Shell script example of how to reproduce simulations using a config file, here 'configs/defaults.toml'
# See documentation for more: <https://pediatriciqphantoms.readthedocs.io/en/latest/notebooks/00_running_simulations.html#Command-Line-Interface-Tool:-make_phantoms>
config='configs/defaults.toml'

make_phantoms $config