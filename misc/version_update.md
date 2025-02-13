To update the version of the package use the folling steps:

* Ad new functionality and check into the repository.
* Change the version number in `Project.toml`.
* Add release notes to `misc/changelog/changelog.md`.
* Push the changes to the repository with the message "Set version to 0.x.y".
* Wait until the bots finish creating the manual etc.
* Then add the comment below to the commit on Githib.

```md
@JuliaRegistrator register

Release notes:

## v0.x.y (February 14, 2025)

- blah
- bla bla
```

