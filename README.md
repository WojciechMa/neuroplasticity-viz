# Neuroplasticity visualisations

Interactive visualisations created as part of a Neuroplasticity course.

## Activity-dependent synaptic refinement

`index.html` simulates how spontaneous retinal waves refine developing connections. A presynaptically gated Hebbian rule potentiates a synapse when its pre- and postsynaptic neurons fire in phase, and depresses and prunes it when they fire out of phase. The page explains every equation in a collapsible section at the bottom.

Live page: https://wojciechma.github.io/neuroplasticity-viz/

## 3D versions (work in progress)

Not linked from the live page yet, and marked `noindex` for search engines. All three run the same learning rule drawn as neurons: in the eye, a retinal ganglion cell sends its axon through the optic nerve, and the axon branches to synapse onto a neuron of each target population.

- `3d/v1/`: action potentials travel along the axons; the thickness of each branch and the size of its synaptic bouton show the synaptic weight.
- `3d/v2/`: v1 plus preset viewpoints, three looks (model colours, Golgi stain, calcium imaging) and an "Explain the parts" mode that describes any part you point at or tap.
- `3d/v3/`: v2 with timing closer to real development and a second eye. Retinal waves are brief bursts that sweep across each retina about once a minute; the clock fast-forwards between waves, and the weights change in steps, once per wave. The two optic nerves meet at the optic chiasm; target A fires with the waves of eye 1 and target B with those of eye 2. Scenarios: waves at different times in the two eyes, synchronised waves, and blocked waves. A timeline shows the firing and the weights over the last five minutes.

They use [Three.js](https://threejs.org/) 0.186.1 (MIT licence), stored in `lib/three/`.

## Run locally

The 2D page is a single file with no dependencies: open `index.html` in a browser. It works offline.

The 3D pages load Three.js as JavaScript modules, which browsers block in files opened from disk. Serve the folder instead, for example with Python:

```
python -m http.server 8000
```

then open http://localhost:8000/3d/v1/, http://localhost:8000/3d/v2/ or http://localhost:8000/3d/v3/.
