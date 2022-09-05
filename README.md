# Covariant amplitudes for decays of meson

The code implements covariant Zemach formalism.
The high-spin expressions and recursive expressions for projectors and orbital tensors are worked out by
Bonn-Gatchina group.

## Groovy implementation

Original author: Misha Mikhasenko

```bash
groovy example_Lj.groovy
groovy example_spin1.groovy
groovy example_spin2.groovy
```

## Python implementation

Original author: Fabian Krinner

The python code is run with
```bash
python makeTensors.py
python makeAmplitude.py
python play.py
```

There are also auxiliary files:
```bash
tensors.np  # Mathematica code for operations with tensors
nonRealtivistic0mp.np # studies of 0-+ decays
```
and 
```
tex2pdf.sh # quick print of the tensor expression in standalone latex.
```

