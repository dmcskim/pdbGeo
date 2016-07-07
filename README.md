# pdbGeo
------
pdbGeo is a structural analysis tool which represents a specified protein region as a directional vector. The angle between the vector and an adjustable plane are then calculated.

This can be used to identify conformational changes among sets of pdb files, like the protein kinase regulatory C-helix conformations. It can also be applied to MD trajectories, with the resulting angles graphed for analysis.

#Python dependencies
------
Please ensure the following Python packages are installed via ``pip``:
- `NumPy <http://www.numpy.org>`
- `MDAnalysis <http://www.mdanalysis.org>`

