## adding packages

in R console:

```R
renv::init()
```

## managing the repo

in a shell:

```
git add -A && git commit -m "init repo"
git push
git tag -a 0.0.1 -m "++" && git push --tags
```