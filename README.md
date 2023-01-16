# xxx

This is the code and configuration for xxx.

You can run this project via [Gitpod](https://gitpod.io) in a web browser by clicking on this badge: [![Gitpod ready-to-code](xxx)](xxx)

* The paper is [here]()
* If you are interested in how we defined our code lists, look in the [codelists folder](./codelists/).
* Developers and epidemiologists interested in the framework should review [the OpenSAFELY documentation](https://docs.opensafely.org)

## project.yaml
The [`project.yaml`](./project.yaml) defines run-order and dependencies for all the analysis scripts. 
**This file should *not* be edited directly**. To make changes to the yaml, edit and run the [`create-project.R`](./create-project.R) script instead.
There is no need to run [`create_project.R`](./create_project.R) if you are simply cloning this repo.

The [`project.yaml`](./project.yaml) is annotated to describe the purpose of each action. 

# About the OpenSAFELY framework

The OpenSAFELY framework is a Trusted Research Environment (TRE) for electronic
health records research in the NHS, with a focus on public accountability and
research quality.

Read more at [OpenSAFELY.org](https://opensafely.org).

# Licences
As standard, research projects have a MIT license. 
