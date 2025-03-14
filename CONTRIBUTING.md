## How to contribute to this project

### Casual contributors

#### Create issues
1. mention the version number of the package
1. write at least one sentence describing the issue
1. share the output of `versioninfo()`
1. indicate if it is a bug report or feature request

### Regular contributors
1. become member of the Telegram channel
1. ask for write access

#### Creating pull requests
1. create a new branch
1. push your changes to the branch
1. create a draft pull request
1. when you think it is ready, mark it as "ready for review" and invite a reviewer
1. address the remarks of the reviewer
1. check that all tests passed
1. merge the pull request, choosing the `squash` option

It is OK to merge your own pull request without review if you are very sure it will not cause any problems. If in doubt, ask on Telegram.

It is OK to push directly to the `main` branch for small changes for example fixing a typo in the documentation.

#### Commit messages
1. shall start with capital letter
1. imperative form

**Examples**
- Fix division-by-zero error
- Improve docstrings
- Implement feature #23
- Add function cool_feature()

#### Coding style
1. use verbs for the names of functions, nouns for types and structs
1. variable and function names lower case, use underscore if needed
1. Types and structs CamelCase
1. indent with four spaces
1. max 120 char per line
1. do NOT use an automated code formatter

#### File formats
1. input files: preferably `.yaml`format
1. CAD geometries: `.obj` files
1. multi-dimensional arrays: `.mat` files
1. CSV files can be used for polars for example, but if the polars are part of wing description better have one `.yaml` file that describes all aspects of the wing

**Open questions** 
- what are `.dat` files?
- shall we use `.bin` files?
- shall we use `.JLD2` files?