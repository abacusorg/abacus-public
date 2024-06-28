## Spack

We would like to make Abacus Spack-installable. This shouldn't be too hard, especially if we switch to Meson.

Currently, we're only using Spack to build a containerized environment for testing (`ci` directory). The same environment can be used to set up Abacus dependencies for production:

```console
$ spack create env abacus ci/spack.yaml
$ spacktivate abacus
```

For self-managed systems, this is a great solution. However, for HPC systems with an existing software stack (e.g. modules), one should probably use that stack instead of Spack. For missing dependencies, one can in theory use Spack, but it doesn't seem easy to get it to use external dependencies instead of its own by default. The closest thing seems to be to install without Spack dependencies:

```console
$ spack install --only=package ...
```

For simple dependencies, this is probably a good solution. More complex dependencies would require doing this for each requirement in the dependency graph, though. In other words, there's no simple way to have Spack only install missing dependencies.

There is a more complex way, and that requires specifying what "external" dependencies are already present (https://spack.readthedocs.io/en/latest/packages_yaml.html#external-packages). Some HPC sites have started distributing `spack.yaml` files that describe these external packages, which can help with this. Whether this is worth it probably depends on the complexity of the missing dependencies.

For the moment, we'll leave some bundled dependencies in `external/`. However, we may want to encourage users to get dependencies via Spack in the future.
