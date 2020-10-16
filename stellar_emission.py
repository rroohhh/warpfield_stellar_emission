import numpy as np
from pathlib import Path
import os
import sys
import re

class StellarPopDataWrapper:
    def __init__(self, data, idx):
        self._data = data
        self._idx = idx

    def __getattr__(self, attr):
        return getattr(self._data, attr)[self._idx]

class StellarPopData:
    def __init__(self, filename):
        with open(filename) as f:
            header = next(f)
            assert header[0] == '#', f"malformed stellar pop file {filename}"

            header_info = header[1:].strip().split(';')
            data = {name: [] for name in header_info}

            for line in f:
                parts = line.split(';')

                for name, part in zip(header_info, parts):
                    if part == "True":
                        part = True
                    elif part == "False":
                        part = False
                    else:
                        part = part.split('\t')
                        part = [float(p) for p in part]

                    data[name].append(part)

            self._data = data

    def __getattr__(self, attr):
        return self._data[attr]

    def cluster_data(self, n):
        idx = np.where(self.id == n)[0][0]

        return StellarPopDataWrapper(self, idx)

class SB99Spectrum:
    def load(filename):
        filename = Path(filename)

        self = SB99Spectrum()
        if (match := re.search(r"([0-9eE\.]+)cluster_(rot|norot)_Z([0-9]{4})_BH([0-9]+).stb99", filename.name)):
            self.data = np.genfromtxt(filename, names=["age", "lambda", "loglum_total", "loglum_stellar", "loglum_nebular"], skip_header=6)

            SB99mass = float(match.group(1))

            if match.group(2) == "rot":
                rotation = True
            elif match.group(2) == "norot":
                rotation = True

            if match.group(3) == "0014":
                Z = 1.0
            elif match.group(3) == "0002":
                Z = 0.15

            BHcutoff = float(match.group(4))

            self.SB99mass = SB99mass
            self.props = (rotation, Z, BHcutoff)
            return self
        else:
            return None


    # dL/dλ in erg / (s * A)
    def population_emission(self, masses, birth_ages, ages):
        ages = np.atleast_1d(ages)
        birth_ages = np.atleast_1d(birth_ages)
        masses = np.atleast_1d(masses)

        all_ages_masses = np.tile(masses, (len(ages), 1))
        all_ages = np.tile(ages, (len(birth_ages), 1)).T - birth_ages
        all_ages_masses[np.where(all_ages < 0)] = 0

        sb99_ages, indices, counts = np.unique(self.age, return_index=True, return_counts=True)
        assert np.all((counts - counts[0]) == 0), "not every age range has the same wavelength grid"
        assert np.all((indices % counts[0]) == 0), "ages were not sorted"

        age_indices = indices[np.argmin(np.abs(np.tile(sb99_ages, (*all_ages.shape, 1)).T - all_ages.T), axis=0).T] // counts[0]

        emission = 10**self.loglum_stellar.reshape((-1, counts[0]))
        emissions = np.sum((emission[age_indices].swapaxes(0,2) * (all_ages_masses / self.SB99mass).T).swapaxes(0,2), axis=1)

        return np.log10(emissions)

    def wavelength_grid(self):
        return self["lambda"][np.where(self.age == self.age[0])]

    def __getattr__(self, attr):
        return self.data[attr]

    def __getitem__(self, attr):
        return self.data[attr]

class StellarEmissionSamples:
    def __init__(self, datadir="./data"):
        self.pop_data = StellarPopData(Path(datadir) / "stellar_pop_data.dat")

        self._load_sb99_spectra(datadir)

    def _load_sb99_spectra(self, datadir: str):
        spectra_dir = Path(datadir) / "sb99_spectra"

        self.spectra = {}
        for spectrum in spectra_dir.iterdir():
            if spectrum.is_file():
                if (spectrum := SB99Spectrum.load(spectrum)) is not None:
                    self.spectra[spectrum.props] = spectrum

    def dump(self, ids = None, ages = None):
        if ids is None:
            ids = self.pop_data.n

        ids = np.atleast_1d(ids)

        wavelength_grid = None
        emissions = []
        ns = []
        done_ages = []

        for n in ids:
            cluster = self.pop_data.cluster_data(n)

            if ages is None:
                cluster_ages = cluster.ages
            else:
                cluster_ages = np.atleast_1d(ages)

            sb99_spectrum = self.spectra[(cluster.rotation, cluster.Z[0], cluster.BHcutoff[0])]
            if wavelength_grid is None:
                wavelength_grid = sb99_spectrum.wavelength_grid()
            else:
                assert np.all(wavelength_grid == sb99_spectrum.wavelength_grid()), "dumping emission of clusters with different wavelength grids not supported yet"

            emissions.append(sb99_spectrum.population_emission(cluster.masses, cluster.birth_ages, cluster_ages))
            ns.append([n for _ in range(len(cluster_ages))])
            done_ages.append(cluster_ages)

        emissions = np.vstack(emissions)
        return (wavelength_grid, np.ravel(np.concatenate(ns)), np.ravel(np.concatenate(done_ages)), emissions)


if __name__ == '__main__':
    emission_samples = StellarEmissionSamples()


    if len(sys.argv) == 1:
        out = sys.stdout
    elif len(sys.argv) == 2:
        out = open(sys.argv[1], "w")

    first = True
    for id in emission_samples.pop_data.id:
        grid, ids, ages, emissions = emission_samples.dump(ids=[id])
        if first:
            first = False
            out.write(f",,{','.join([str(g) for g in grid])}\n")

        for i, a, e in zip(ids, ages, emissions):
            out.write(f"{i},{a},{','.join([str(v) for v in e])}\n")

    if out != sys.stdout:
        out.close()
