# Simulation of the Evolution of Primordial Fluctuations into the Cosmic Web

**Author:** María Constanza Lepe  

This Python project simulates the gravitational evolution of primordial density fluctuations — tiny variations in matter density that originated during the early universe (inflation) — into the large-scale structure we observe today: galaxy filaments, clusters, and cosmic voids.

The code:
- Generates a **3D Gaussian random field** following a power spectrum \( P(k) \propto k^{-3} \).
- Solves the **Poisson equation** in Fourier space to compute the gravitational potential.
- Uses a **Leapfrog integration scheme** to simulate the growth of structures.
- Incorporates a simple **cosmic expansion** term for a more realistic evolution.
- Produces visualizations of density fields and power spectra.

> ⚠️ This is a **toy model** intended for educational and portfolio purposes, not a fully accurate cosmological simulation.

---

## 🔬 Scientific Background

In cosmology, large-scale structure forms from initial fluctuations that grow under gravity.  
This simulation starts from a Gaussian random field in the early universe and shows how overdense regions collapse, forming structures similar to the **cosmic web** seen in galaxy surveys such as the **Sloan Digital Sky Survey (SDSS)**.

---

## 🚀 How to Run

### 1. Install dependencies
```bash
pip install numpy matplotlib


2. Run the script

python cosmic_web.py

You will see:

-Initial density field slice.
-Power spectrum in Fourier space.
-Gravitational potential map.
-Final evolved density after toy gravitational dynamics.



📄 Requirements

Python 3.9+

NumPy

Matplotlib



📌 Notes

Grid size can be changed in grid_size variable, but higher values require more RAM.

For professional cosmological simulations, methods like N-body codes or Zel’dovich approximation are recommended.



📚 References

Peebles, P. J. E. Principles of Physical Cosmology

Springel, V. et al., The Millennium Simulation

SDSS SkyServer


