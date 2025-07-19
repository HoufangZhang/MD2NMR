from setuptools import setup, find_packages

setup(
    name="MD2NMR",
    version="0.0.5",
    author="Tiejun Wei",
    description="Tools for calculating NMR relaxation observables (R1/R2/NOE/T1/T2/Tau_c) directly from MD trajectories. Initially written for calculations regarding nucleosome simulations but can be extended for other proteins/complexes. version 0.0.5: 1.Updated with downloading testing MD traj on zenodo: https://zenodo.org/record/7806382 and a new readme file for testing different MD traj. 2. Updated a new function allowing user to select residues with use_chain_id option in config.py file. ",
    url="https://github.com/DerienFe/MD2NMR",
    packages=find_packages(),
    install_requires=["numpy == 1.23.4",
                      "pandas == 1.5.1",
                      "scikit-learn == 1.1.3",
                      "scipy == 1.9.3",
                      "mdtraj ==1.9.7",
                      'setuptools == 65.5.0'],
)