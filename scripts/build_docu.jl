# build and display the html documentation locally
# you must have installed the package LiveServer in your global environment

using TestEnv; TestEnv.activate()
using LiveServer; servedocs(launch_browser=true)
