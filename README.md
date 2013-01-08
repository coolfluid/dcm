**Coolfluid3 plugin "sdm":**

Spectral Difference Method for solving Partial Differential Equations

**Installation:**

  + Create a plugins directory, and clone the sdm sources inside:

```
mkdir -p $CF3_PLUGINS_DIR
cd $CF3_PLUGINS_DIR
git clone https://github.com/coolfluid/sdm.git $CF3_PLUGINS_DIR/sdm
```

  + Rerun cmake in the coolfluid3 build directory:

```
cd $CF3_BUILD_DIR
cmake .  -DCF3_PLUGINS_DIR=$CF3_PLUGINS_DIR -DCF3_PLUGIN_SDM=ON
```
