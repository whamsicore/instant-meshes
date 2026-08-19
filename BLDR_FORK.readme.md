<!-- BLDR_FORK.readme.md -->

# Instant Meshes Fork (BLDR)

Vendored fork of [wjakob/instant-meshes](https://github.com/wjakob/instant-meshes) for Dharmaland Retopo Lab.

**Name**: Instant Meshes Fork — so it is obvious this is our controllable fork, not a mysterious “Instant Meshes” product mode.

## Build

```bash
cd vendors/integrated/instant_meshes_fork  # after vendor_sync into BLDR-CLUB
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_POLICY_VERSION_MINIMUM=3.5
cmake --build build -j "$(sysctl -n hw.ncpu)"
```

CLI: `build/instant_meshes_fork.app/Contents/MacOS/instant_meshes_fork`

## Field export (BLDR)

```bash
# Orientation solve only — write Q/N/V CSV (no remesh extract)
instant_meshes_fork --field-only --export-field field.csv sealed.obj

# Remesh + also dump field
instant_meshes_fork -o out.obj -f 200000 --export-field field.csv sealed.obj
```

Stroke GLB for Forge View:

```bash
Blender --background --python scripts/dharmaland/constrained_retopo/export_im_orientation_field_strokes.py -- \
  --field-csv field.csv --base-obj sealed.obj --output strokes.glb
```

## Rules

1. Never publish remesh to catalog unless topology hole/manifold audit passes.
2. No closed voxel remesh for character retopo.
3. Fail the stage if seal prep fails — do not ship cheese-grater meshes.
