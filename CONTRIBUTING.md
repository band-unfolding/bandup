# Contributing to BandUP
#### BandUP: Band Unfolding code for Plane-wave based calculations             
###### Copyright (C) 2013-2017 Paulo V. C. Medeiros - pvm20@cam.ac.uk 
###### Please visit <http://www.ifm.liu.se/theomod/compphys/band-unfolding>

<!-- =========================================================== -->
## Introduction          

**First and foremost, thank you for your interest in BandUP and for
considering to contribute to it**. We greatly appreciate the valuable
feedback given by all of those who have contacted us to suggest
features, provide benchmarks for our tests, report issues, etc., and
we welcome collaboration from all those who would like to do some
coding as well.

This document will guide you through the the workflow and conventions
adopted in the development of BandUP. These have been established to
promote consistency and minimise the occurrence of conflicting and/or 
ill-behaving code, thus making collaboration easier for everyone.

The following instructions assume you have basic knowledge of Git.
If you don't know any Git but wish to contribute anyway, then please
get in touch with me by email. Contributions are always welcome.

<!-- =========================================================== -->
## A. Workflow options

To collaborate with the BandUP project, **we recommend** that you 
create an account on GitHub if you don't already have one. This will
allow you to adopt the **Fork-and-Branch workflow**, which we prefer.
If using a GitHub account is not a viable option for you, then you may
also adopt the *Feature Branch workflow*, although this is *not the 
approach we'd normally recommend for this project*. We discuss both
workflows, within the context of BandUP, in sections A.1 and A.2. For
more general information about git, forking repos, and workflows,
please visit, eg., 
<https://help.github.com/articles/fork-a-repo> and 
<https://www.atlassian.com/git/tutorials/comparing-workflows>.

Regardless of the picked workflow, however, BandUP adopts the following
**"dedicated branch convention"**:
<p align="center">
          **Development takes place in dedicated branches.**
</p>
More specifically, **the official repo's `master` and `devel` shall
only be modified by BandUP's main developers**.

  * `master` shall only be used for releases, and shall receive merges
    from `devel` only  
  * Other **branches shall be created by branching off from `devel`**,
    and may be merged back into `devel` once ready

Rationale: Most users will only ever download `master`, and they 
naturally expect it to contain usable, trustworthy code. It is thus
safer to have changes made in dedicated branches. BandUP's main 
developers will then revise the new/modified code and, after approval, 
merge such branches into `devel` (for integration with other branches 
that may have been developed parallelly). Code in `devel` will be 
merged into `master` once a new release is ready.

The following instructions explicitly enforce the dedicated branch
convention. When following them, beware that **any changes you may
have made to the tracked files on `devel` will be overwritten**.

### A.1: Fork-and-Branch workflow (preferred method)

#### First-time setup

1. Fork your own copy of the official BandUP repo on GitHub
2. Make a local clone of *your fork*. You will work on this clone for
   the rest of the workflow
3. Define a new remote `upstream` pointing to the official repo:

        git remote add upstream https://github.com/band-unfolding/bandup.git

#### Syncing your local `devel` with `upstream`

Get your local `devel` up-to-date with the the official repo: 

                git fetch upstream
                git checkout devel
                git reset --hard upstream/devel

**NB.:** Any changes you may have made to the tracked files in 
**your local `devel` will be overwritten**.

#### Working on a dedicated branch        
        
  * If you don't already have a branch designated for what you
    plan to do, then create one. **Always branch off from `devel`**:
                  
                git checkout -b your_branch devel

  * If you already have a branch, then make sure to be using it:
      
                git checkout your_branch
              
    and then get it up-do-date with `devel`:

                git rebase devel
 
Now you can safely modify the code.

#### While modifying/adding code
   
  * **Stick to the coding conventions** adopted in the file you are
    editing. If adding a new file, use one of the existing files as
    template. Aim at a **maximum linewidth of 80 characters**, so
    that code can easily be visualised in small screens.
  * Don't be shy to **commit your changes whenever reasonable**.
    Smaller commits make it easier to traceback and fix eventual bugs.
  * **Use meaningfull commit messages**. Other people will need to
    understand it!

#### Once you are done
 
   1. Make sure your code can be compiled and has no obvious bugs.
   2. Commit your final changes.
   3. Rebase the remote `devel` from the official repo into your branch
     (new code may have been pushed while you were working):
     
                git fetch upstream
                git checkout your_branch
                git rebase upstream/devel
      
   4. Make sure again that the resulting code can be compiled and has
      no obvious bugs.
   5. Push your changes *to your remote fork*:

                git push -u origin your_branch
                
   6. Create a pull request from your branch into the official repo's
      `devel`. You do this via GitHub's online interface (the process
      is quite self-explanatory).
      

### A.2: Feature Branch workflow (less preferred method)
#### First-time setup

If you haven't already done so, you must make a local *clone* of 
BandUP's official repo (and not just download it). You can do so by
typing (no password required):
   
            git clone https://github.com/band-unfolding/bandup.git
                
#### Syncing your local `devel` with `origin`
If working on a pre-existing clone, get your local `devel` up-to-date
with the remote one:

                git fetch origin
                git checkout devel
                git reset --hard origin/devel

**NB.:** Any changes you may have made to the tracked files in 
**your local `devel` will be overwritten**.
        
#### Working on a dedicated branch        
        
  Follow the same instructions provided for this task in the case of
  the Fork-and-Branch workflow. After this, you can safely modify the
  code.

#### While modifying/adding code
   
  Follow the same instructions provided for this in the case of the
  Fork-and-Branch workflow.

#### Once you are done
 
   1. Make sure your your code can be compiled and has no obvious
      bugs.
      
   2. Commit your final changes.
      
   3. Rebase the remote `devel` into your branch (new code may have 
      been pushed to the remote while you were working on your local
      repo):

                git fetch origin
                git checkout your_branch
                git rebase origin/devel
      
   4. Make sure again that the resulting code can be compiled and has
      no obvious bugs.
      
   5. **Email us** asking for membership to the BandUP organisation
      if you haven't yet done so. This will enable you to push your
      branch to our GitHub remote. *This step is one of the reasons 
      why we strongly recommend the Fork-and-Branch workflow*.
      
   6. Finally, push your branch:

                git push -u origin your_branch

<!-- =========================================================== -->
## B: After you submit your changes 
                
After BandUP's main developers are aware of your new code, they will check your proposed changes/additions. If there are any issues, you 
will be contacted. If they are approved, your code will be merged into
the official repo's `devel` for further integration with any other
eventual existing development. All new features, fixes, etc. will be
merged from `devel` into `master` when the new code is tested, 
approved, and a new release can thus be created.               
                
<!-- =========================================================== -->
## C: For BandUP's main developers 

As of BandUP v3.0.0-beta.1, the idea is to make sure that commits to
`master` are tagged according to ***semantic versioning***, and that 
the tags used are consistent with the code version displayed at
runtime. To this end, the following workflow must be adopted to push
new code to `origin/master`:      
      
   1. As previously discussed, `master` shall only receive merges from 
      `devel`. Changes from other branches must thus be integrated into
      `devel`. Before merging into `devel`, please certify that the
      code to be merged:
      
      * Can be compiled
      * Has been tested
      * Has no obvious bugs.

   2. Merge into `devel`, one by one, the user branches/pull requests
      that are ready for such.
      
      * **Do not rebase**.
      * Use **`--no-ff`**
      * Use your best judgement to solve conflicts. If you've got
        permission to push to `origin/master` and `origin/devel`, that
        means I trust you on this.
      * After each merge, make sure that the resulting code can be
        compiled and has no obvious bugs.
   
   3. Before making the final commit to `devel` (i.e., the commit to be
      merged into `master`), update the value of `tag_for_push` 
      (`constants_and_types_mod.f90` file, under `src`). 
      Follow the **semantic versioning 2.0** set of conventions. 
      See <http://semver.org>.
           
   4. Now you can commit any new changes (such as the one made 
      in the previous step) to your local `devel`.
           
   5. Merge `devel` into `master`. Use **`--no-ff`**.
           
   6. **Tag the merge commit made in the previous step**. The created
      tag must be **(a)** identical to the value attributed to
      `tag_for_push` (see item 3 above), and **(b)** of the
      **annotated** type. See 
      <https://git-scm.com/book/en/v2/Git-Basics-Tagging>.
      * I use a script for (5) and (6). It retrieves the value of 
        `tag_for_push` from the file mentioned in (3), merges `devel`
        into `master`, and then creates the appropriate annotated tag
        with a message I choose.
           
   7. **Double-check everything**. This is the last chance to 
      correct mistakes in a straightforward manner.
           
   8. Finally, push `master`, `devel`, and the new tag to the remote.

<!-- =========================================================== -->
## Closing remarks         

Please feel free to get in touch with us if you need assistance with
anything that has been discussed here. We also encourage you to come
forward with suggestions of new features, workflow and code improvements, corrections to the information presented here, etc. Any
contribution is welcome. 

Once more, thanks for your interest in BandUP and for considering to
contribute to it.
            
**Have fun!**  
Cheers,  
Paulo.
