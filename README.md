Computational model of Base-Excision-Repair (BER) mechanism
===========================================================

This repository contains the implementation of the computational model of
BER mechanism, as described in [paper]


Requirements
------------

* KaSim 3.5
* Python 2.7
* simplejson


Usage
-----

Main entry point is the "run.py" file. It requires the user to specify
the run mode as it's first command line argument, e.g:

```
$ python run.py [mode]
```

where `mode` is one of the following:

* `synth`, which synthesizes sample DNA strands. It supports optional 
named argument `-descriptor` that points to the JSON settings file.
* `simulate`, which simulates the BER mechanism.  It supports optional 
named argument `-descriptor` that points to the JSON settings file.
* `experiment`, which runs a series of experiments. It requires an 
additional argument that points to a JSON file that describes the 
experiment.

All values defined in the JSON files can be overridden in command 
line.


`synth` and `simulate` mode have a set of predefined default values. 
The supplied JSON file (`-descriptor`) needs not to define every value, 
instead, missing values are filled from the default dictionary.

JSON files are separated into `"meta"` and `"data"` top-level keys. 
`"meta"` defines the "meta" values, used to direct output, etc. "data"
is instead used to define values used directly in the simulation.

`synth`
-------

* `meta:temp_dir` 

	defines the temporary directory used in the synthesis

* `meta:out_dir` 

	defines the output directory

* `meta:out` 

	specifies the filename of the synthesized DNA

* `data:bp_per_strand` 

	defines the number of base pairs per strand

* `data:strands` 

	defines the number of strands in the output

* `data:weights` 

	defines the dictionary of weights used when generating the DNA

`simulate`
----------

* `meta:temp_dir` 

    defines the temporary directory used in the simulation

* `meta:out_dir` 

    defines the output directory

* `meta:out` 

    specifies the filename of the trace of the simulation

* `meta:dna_file` 

    specifies the filename of the DNA file, on which the BER is performed

* `meta:time` / `meta:events` 

    specifies the amount of "simulation time"/events used for the 
    simulation

* `meta:plot_points`

    specifies the number of points in the trace of the simulation

* `meta:model`

    is an array that defines the components used in the simulation:

    * `damage` -- damages the DNA file after 1000 events
    * `sliding` -- sliding action for the DG actor
    * `DG`, `APE1`, `POLb`, `LIG3`, `XRCC1`, `PNKP` -- actions for the
    respective actor
    * `xrcc_dimer` XRCC1 dimerization action

* `meta:sanity`

	used for debugging. Outputs weakly compressed stories of some 
	observables.

* `data:...` 

	specifies the kinetic rates used in the simulation -- see 
	rates.xlsx

`experiment`
------------

`meta` fields are joined with the ones defined before (in `simulate`).

* `meta:experiment_name` 

	defines the name of the experiment series

* `meta:runs` 
	
	specifies the number of times the rates are calculated and 
	simulations are performed

* `meta:parallel`

	specifies the number of simulations to perform in parallel

* `data:constants` 
	
	specifies the dictionary of constants used in calculation of
	rates

* `data:fields`

	specifies the dictionary of fields to generate the rates for. The
	keys of the dictionary should match the names of the rate 
	variables. The value of each key should be another dictionary with
	two components `type` and `value`.

	Currently supported types are: `constant`, `random_uniform`, 
	`Kd_constraint`.

	* `constant` copies the value of the `value` key as the value of 
	the rate variable
	* `random_uniform` generates a uniform random value in the range
	specified in the `value` key. The range is defined as a JSON array
	of three values, the first one is either `"int"` or `"float"`, 
	followed by the lower bound, and the upper bound.
	If the first value is missing, `"float"` is used. You can also use
	the name of a variable defined in the `constants` section, if the
	variable defines an array of the form specified before.
	* `Kd_constraint`, used to impose the K_d constraint on variables
	that are generated. The value field should be an array, with the
	rate variable name used as k_off in the first position, followed
	by Kd ratio.
