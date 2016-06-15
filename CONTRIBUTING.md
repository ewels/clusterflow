# Cluster Flow: How to Contribute

### Making an Issue
First - before you start working on a change to the repository, please make
sure that there are no exisiting
[issues](https://github.com/ewels/clusterflow/issues) relating to
whatever change you intend to make. If there aren't, please create one
so that others know that you're working on something.

### Workflow
The workflow for adding to this repository should be as follows:

1. [Create an issue](https://github.com/ewels/clusterflow/issues)
   describing what you intend to work on
2. Fork the [development branch](https://github.com/ewels/clusterflow/branches) of the repository to your own GitHub account
	1. Any changes you make will be merged into this fork. If you make changes to the stable master branch instead, this will be a much more painful process.
3. Make your changes. Remember to note these in the `README.md` changelog.
5. Submit a Pull Request describing your changes. I will review your code and merge.

### Retrospective Workflow
This is all well and good if you haven't already started hacking the code. If you downloaded a static version of the code and made your changes, this is the ideal workflow:

1. Register with [github.com](https://github.com/) if you haven't already
	1. There are excellent github tutorials about [forking repositories](https://help.github.com/articles/fork-a-repo/) and creating [pull-requests](https://help.github.com/articles/using-pull-requests/) .
2. Fork the [development branch](https://github.com/ewels/clusterflow/branches) of Cluster Flow to your own GitHub account
3. Pull this fork to your system using [`git clone`](https://help.github.com/articles/fetching-a-remote/)
4. Replace the downloaded files with your modified versions.
5. Push your updates using `git commit -a -m "Message"` and `git push`.
6. Create a [pull request](https://github.com/ewels/clusterflow/pulls) so that the changes can be merged with the development branch.


### Getting Help
If you have any queries, please get in touch with [Phil Ewels](https://github.com/ewels). Thanks for contributing!
