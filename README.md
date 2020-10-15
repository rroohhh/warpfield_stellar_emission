# Warpfield stellar emission
This code generates examples of the stellar emission from warpfield models.
Each example can contain multiple generation of clusters at various ages.

## Usage
To simply dump all emission data for all models at every age deemed "interesting" by warpfield, run
```bash
python3 stellar_emission.py output_filename
```
This will generate a table, where the first column contains a model id, the second column the age of the oldest cluster at the point the emission was calculated and all further columns the luminosity `log10(dL/dλ)` in `erg/(s A)`.
The first row contains the wavelength in Ångström for each of the columns containing the luminosity.

The table at `data/stellar_pop_table.csv` lists the birth time and mass for each of the ids.

### API
To work with the provided population examples, use the `StellarEmissionSamples` class.

Using the `dump` function you can generate the emission for specific models at specific ages:
```python
from stellar_emission import StellarEmissionSamples

em_samples = StellarEmissionSamples()

wavelength_grid, model_ids, ages, emission_data = em_samples.dump(ids=[1,3,100], ages=[1e6, 10e6, 30e6])
```
In `emission_data` this will return the emission data of the models 1, 3 and 100 each at a age of 1e6, 10e6 and 30e6.

For complete control over the generated population you can use the `SB99Spectrum` class:
```python
import numpy as np
from stellar_emission import SB99Spectrum

s = SB99Spectrum.load('data/sb99_spectra/1e6cluster_rot_Z0014_BH120.stb99')

wavelength_grid = s.wavelength_grid()

emission_data = s.population_emission(masses=[1e5, 1e6], birth_ages=[0, 10e6], ages=np.linspace(0, 30e6))
```
Here `emission_data` will contain the combined emission of two clusters, one with a mass of 1e5 and one with a mass of 1e6, with the first created at time 0 and the second created at 10 Myr at different observation times specified by `ages`.

To explore the provided populations, you can use the `StellarPopData` class:
```python
from stellar_emission import StellarPopData

pop_data = StellarPopData('data/stellar_pop_data.dat')

ids = pop_data.id

cluster_data = pop_data.cluster_data(ids[0])

masses = cluster_data.masses
birth_ages = cluster_data.birth_ages
Z = cluster_data.Z
```

Here `ids` is a list of the available model ids, `masses` contains a list of the masses of the star clusters created at `birth_ages` over the time evolution of the model.
