# Barebones environment that can be used on a cluster
It is much like the main env where it has all the code from this repo `dev`ed but it does not include GLMakie or MathLink both of which rely on non-julia outside code that cannot be easily installed on a cluster.
This env is intended for longer/bigger runs that will be ran on the cluster, subdirs in `runs` corresponding to one run each including some of their code for later reference and their results should be written there too but not included in the repo as they will be large files.
