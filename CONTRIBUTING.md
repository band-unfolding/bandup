######################################################################
## Contributing to BandUP            
######################################################################

First and foremost, thank you for your interest in BandUP and for
considering to contribute to it.

<!-- =========================================================== -->
### Semantic versioning and other conventions

The idea here (as of v3.0.0-beta.1) is to automate semantic 
versioning -- fully for any development branch, and at least
partially for master. To this end, the following workflow must 
be adopted (assumes you know Git):
    
   1. No modification in the code should be made in the master 
      branch. Use the "devel" branch instead, or create another 
      branch if appropriate. If working on an existing branch, 
      then, before doing anything, make sure to pull in any 
      eventual new changes from master (which should be up-to-
      date with the remote master). Once you are happy with the
      new/modified code and convinced that it has no obvious
      bugs, then go to (2).
      
   2. Make sure again the code to be merged has been tested and
      has no obvious bugs.
      
   3. Before making the final commit to the branch, the value  
      attributed to the constant "tag\_for\_push" (in the file 
      "constants_and_types_mod.f90" under "src") must be
      updated. This must be done in accordance with the
      semantic versioning 2.0 set of conventions. 
      See <http://semver.org>.
           
   4. After making sure to comply with the items above, commit 
      any new changes -- such as the one you did in (3).
           
   5. If you don't have permissions to push to master, then 
      push your branch, send me a notification, and stop here.
      If you *can* push to master, then merge your branch into
      master and continue. Use "--no-ff" unless you have only
      one commit to merge.
           
   6. Before pushing the merge commit, it must be tagged.
      Moreover: The created tag must be **(a)** identical to the 
      value attributed to "tag\_for\_push" in (3), and **(b)**
      an annotated tag.
      See <https://git-scm.com/book/en/v2/Git-Basics-Tagging>.
      * I use a script for (5) and (6). It retrieves the value of 
        "tag_for_push" from the file mentioned in (3), merges
        the changes made in the branch into master, and then
        creates the appropriate annotated tag with a message I
        choose.
           
   7. Now, double-check everything. This is the last chance to 
      correct mistakes in a straightforward manner.
           
   8. Finally, push to the remote repo.
        
By doing as described above, you will help me greatly and the 
semi-automated semantic versioning scheme should work. 
    
If you don't know Git, but have implemented something you think
would help improve BandUP and you'd like to share it with the
community, then please get in touch with me by email. But I
do recommend you learn Git if you are modifying BandUP source!

##### Have fun!
Cheers,  
Paulo.
