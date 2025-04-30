# Welcome to HCBC Platform

Welcome to our Platform guidelines web-page. Here analysts and developers will find guidelines on how to work with our most common environments.

We are developing platforms for each analysis type we have experience with at HCBC.

Most analyses will follow the similar trajectory for set-up. We will note where the trajectories diverge with tables where you will select the appropriate analysis.

## Set up the package

<details>
  <summary>O2 instructions- Click to expand!</summary>

Log onto O2 via the command line and check two things (first-time only): 

   * Remove `bcbio` from you `PATH` by commenting the line in your `.bashrc` if you have it
   * Remove any path you load using the R env variables that could be in your `.Rprofile` or `.bashrc`

Go to the [O2 Portal](https://o2portal.rc.hms.harvard.edu/) and select `HMS-RC Application`, then `RStudio Environment`
<p align="center"><img src="./img/O2_portal_page_with_label.png" width="700"></p>

Start Rstudio with using your desired partition, memory, core and time directives. Add these modules to "Modules to be loaded":
  ```
  git/2.9.5  gcc/9.2.0 imageMagick/7.1.0 geos/3.10.2 cmake/3.22.2 R/4.3.1 fftw/3.3.10 gdal/3.1.4 udunits/2.2.28  boost/1.75.0
  ```
   Click "Launch"

<p align="center"><img src="img/O2_portal_request.gif" width="700"></p>

Open RStudio by clicking on the "Connect to RStudio Server"

<p align="center"><img src="img/O2_connect_to_RStudio_with_label.png" width="700"></p>

When the session is started, set your library path by typing this command in your console Rstudio window in order to be able to load `bcbioR`:

```
.libPaths("/n/app/bcbio/R4.3.1")
```
</details>

</p>
Next, load `bcbioR` with:

```
library(bcbioR)
```

Check the package version of `bcbioR` using:

```
packageVersion("bcbioR")
```

Make sure the version is **0.4.*** or later.

<p align="center"><img src="./img/Version_check_bcbioR.gif" width="700"></p>

> Note: If you are working in your local environment, install `bcbioR` with
> ```
> devtools::install_github("bcbio/bcbioR",build_manual = TRUE, build_vignettes = TRUE)
> ``` 

## General Project

Use the app in the dropdown below to use the HBC app to set up a project's name. This name will be used for O2, FAS, GitHub and Dropbox. This name is likely already defined in the Trello card.

<details>
<summary><b>Click here to use the HBC app for naming a project</b></summary>
<p align="center"><iframe src="https://hcbc.connect.hms.harvard.edu/content/8cd62872-0ec9-4905-8920-c745d2375758/" width="800px" height="1000px" data-external="1"></iframe></p>
</details>

This set up needs [bcbioR](https://github.com/bcbio/bcbioR) and [usethis](https://usethis.r-lib.org) packages. If you are working on the O2 Portal `usethis` is within `/n/app/bcbio/R4.3.1`, so you will not need to install it. Also, `usethis` is a dependency of `bcbioR`, so if you are working locally it should have come along with the download of `bcbioR`.

### Create the Rstudio project

Assign the path that you will be using as the path for your project to the object `project_path` into the `/n/data1/cores/bcbio/PIs` space on O2:

```
project_path <- "/n/data1/cores/bcbio/PIs/PI_name/lastname_postdoc_rnaseq_human_heart_hbc00000"
```

Now, we will need to create a directory using this path to put our analysis in. If this directory already exists, then you can skip this step. The directory creation can be done with this command:

```
dir.create(project_path)
```

Now, we can open an Rproject for our analysis in this path using:

```
usethis::proj_activate(project_path)
```

> Note: This will restart the session in the project directory. This restart will clear the `.libPaths("/n/app/bcbio/R4.3.1")` and `library(bcbioR)` that we used earlier, so we will need to re-do them in the following steps. 

## Using the template reports

Many analyses have template reports that you can use. You can use these by using the approriate `bcbioR::bcbio_templates()` command from the table below:

| Type of Analysis | `bcbioR::bcbio_templates()` command |
|:---:|:---|
| Bulk RNA-seq ![](https://img.shields.io/badge/status-stable-blue)| `bcbioR::bcbio_templates(type="rnaseq", outpath="reports")` |
| Single-cell RNA-seq ![](https://img.shields.io/badge/status-beta-yellow) | `bcbioR::bcbio_templates(type="singlecell", outpath="reports")` |
| ChIP-Seq ![](https://img.shields.io/badge/status-beta-yellow) | `bcbioR::bcbio_templates(type="chipseq", outpath="reports")` |
| CellChat ![](https://img.shields.io/badge/status-draft-grey)| **Under development:** `bcbioR::bcbio_templates(type="singlecell_delux", outpath="reports")` |
| COSMX ![](https://img.shields.io/badge/status-draft-grey)| `bcbioR::bcbio_templates(type="spatial", outpath="reports")` |
| DNA Methylation | **Under development** |

## Setting up your workspace in O2

We will now add the `.libPath()` that is appropriate for our type of analysis. You can use the table below to determine which `.libPath()` is appropriate for your analysis:

| Type of Analysis | `.libPath()` command |
|:---:|:---|
| Bulk RNA-seq| `.libPaths("/n/app/bcbio/R4.3.1_rnaseq")` |
| Single-cell RNA-seq | `.libPaths("/n/app/bcbio/R4.3.1_singlecell")` |
| ChIP-Seq | `.libPaths("/n/app/bcbio/R4.3.1_chipseq")` |
| CellChat | `.libPaths("/n/app/bcbio/R4.3.1_cellchat")` |
| DNA Methylation | `.libPaths("/n/app/bcbio/R4.3.1_methylation")` |


Since our loaded libraries were wiped when we created a new project, we will need to reload `bcbioR`:

```
library(bcbioR)
```

Once again, let's just double check our version of `bcbioR` to make sure we are using version **0.4.** or later.

```
packageVersion("bcbioR")
```

Next, we need to set the library that `bcbioR` will use each time you open up this RStudio project. By using the `bcbioR::use_library()`, you will be adding the path to your project `.Rprofile` and this will be sourced each time when you open this RStudio project.

| Type of Analysis | `bcbioR::use_library()` command |
|:---:|:---|
| Bulk RNA-seq| `bcbioR::use_library("/n/app/bcbio/R4.3.1_rnaseq")` |
| Single-cell RNA-seq | `bcbioR::use_library("/n/app/bcbio/R4.3.1_singlecell")` |
| ChIP-Seq | `bcbioR::use_library("/n/app/bcbio/R4.3.1_chipseq")` |
| CellChat | `bcbioR::use_library("/n/app/bcbio/R4.3.1_cellchat")` |
| DNA Methylation | `bcbioR::use_library("/n/app/bcbio/R4.3.1_methylation")` |

Now, we will use `bcbioR` to set-up the directory structure that we will be using for our analysis using the following command:

```
bcbioR::bcbio_templates(type="base", outpath=".", org="hcbc")
```

## Setting up GitHub and RStudio

Now, we will connect O2 with GitHub. First, check in your Home directory if a `.gitconfig` file exists. ***You should only need to do this once.*** The contents should look like:

```
[credential]
	helper = store
[user]
	email = your_email@gmail.com
	name = your_GitHub_username
```

Where `your_email@gmail.com` is the e-mail associated with your GitHub account and `your_GitHub_username` is your GitHub username.

If you do not have a `.gitconfig` file in your Home directory, follow the instructions in the dropdown below to create a `.gitconfig` file. 

<details>
<summary><b>Click here for instructions for making a <code>.gitconfig</code> file</b></summary>
<br>In order to make a <code>.gitconfig</code>, you will need to run the following command:<br><br>
<pre>
usethis::use_git_config(user.name = "your_GitHub_username", user.email = "your_email@gmail.com")
</pre>
You will replace <code>your_GitHub_username</code> with your GitHub username and <code>your_email@gmail.com</code> with the e-mail associated with your GitHub account.<br><br>
<blockquote>
Reminder: If you already have your GitHub username and e-mail associated with your GitHub account in your `.gitconfig` file then you can skip this step.
</blockquote>
Click on circular arrow in the top right of the <code>Files</code> tab called <code>Refresh file listing</code>. You can check that you have successfully made this <code>.gitconfig</code> file by looking for the hidden file in your home directory. The content should look like:<br><br>
<pre>
[credential]
	helper = store
[user]
	email = your_email@gmail.com
	name = your_GitHub_username
</pre>
<blockquote>
Note: In order to see hidden file in your file browser on the O2 Portal, you will need to click on the cog in the top right of the <code>Files</code> tab and select <code>Show Hidden Files</code>.
</blockquote>
<hr />
</details>

### Getting the Git tab

Now, we would like to get the Git tab into our Workspace Browser (where `Environment`, `History`, `Connections` and `Tutorial` tabs are located). We show this transition below:

<p align="center"><img src="./img/Git_pane_RStudio.png" width="1000"></p>


In order to do this we will use the following command:

```
usethis::use_git()
```

This command is the command one would usually use in order to make a commit. However, if we choose to not make a commit and then restart R, it will give us our Git tab in the Workspace Browser. This command should return:

```
✔ Setting active project to "/home/wig051/git_demo".
✔ Initialising Git repo.
✔ Adding ".Rdata", ".httr-oauth", ".DS_Store", and ".quarto" to .gitignore.
ℹ There are 6 uncommitted files:
• .gitignore
• .Rprofile
• DataManagement-Checklist.pdf
• information.R
• README.md
• reports/
! Is it ok to commit them?

1: Yes
2: Negative
3: Not now
```

We do not want to commit these files yet, so select one of the answers like `Negative`, `No way` or `No`.

> Note: It is unclear why Git gives us three choices. The possible choices are re-arranged each time and have different text. So your option may vary a bit from these, select one of the appropriate options for **NOT** commiting.

Next, you will be prompted as to whether you want to restart R:

```
A restart of RStudio is required to activate the Git pane.
Restart now?

1: For sure
2: Negative
3: Not now
```

We will need to restart R in order to get the Git tab in our R Studio, so select `For sure`, `Yeah` or some other option for answering in the affirmative. 

### Creating the first commit

Now, we are going to create our first commit. In order to do this, we need to:

1. Navigate to the Git tab in our Workspace Browser
2. Left-click to check the "Staged" checkboxes next to `README.md` and `.gitignore` to select the files that we would like to add to our first commit
3. Left-click <kbd>Commit</kbd> in the Git tab
4. Add a message for our commit
5. Left-click <kbd>Commit</kbd> underneath the textbox where we added the message for our commit
6. Close both Git windows

These steps are summarized in the GIF below:

<p align="center"><img src="./img/Initial_commit.gif" width="1000"></p>

### Pushing our initial commit

Now we will use the function to push these changes to GitHub with the following command:

```
usethis::use_github(org="hbc", private=TRUE)
```

If you do not have a valid GitHub token, this will give you an error (see [below](#expired-or-non-existent-github-token)). If you already have a valid GitHub token, you will get this text in your console and it will push this repository to the HBC GitHub with the following output:

```
ℹ Defaulting to "https" Git protocol.
✔ Setting active project to "/home/wig051/git_demo".
✔ Creating private GitHub repository "hbc/git_demo".
✔ Setting remote "origin" to "https://github.com/hbc/git_demo.git".
✔ Pushing "master" branch to GitHub and setting "origin/master" as upstream branch.
✔ Opening URL <https://github.com/hbc/git_demo>.
```

If the push is successful, then it will look like this GIF below:

<p align="center"><img src="./img/GitHub_initial_push.gif" width="1000"></p>

> Note: You might get a GitHub 404 error page (see image below) when you do your first push to GitHub. Just refresh the page in your browser and it should be resolve itself.
> <p align="center"><img src="./img/GitHub_404_error_with_label.png" width="700"></p>

#### Expired or non-existent GitHub token

However, if your token is expired or this is your first time using GitHub from O2, then you will get this message:

```
ℹ Defaulting to "https" Git protocol.
✔ Setting active project to "/home/wig051/git_demo".
Error in `usethis::use_github()`:
✖ Unable to discover a GitHub personal access token.
ℹ A token is required in order to create and push to a new repo.
☐ Call usethis::gh_token_help() for help configuring a token.
Run `rlang::last_trace()` to see where the error occurred.
```

If this is the case, then you will need to create a new GitHub token, click on the dropdown box below for instructions on how to create a new GitHub token.

<details>
<summary><b>Click here to see how to create a GitHub Token</b></summary>
<br>The first thing we will need to do when setting up our GitHub token is to set-up our credential helper. We can do this by left-clicking on the <code>Terminal</code> table in the Console Window. You may need to load Git in the Terminal using:<br><br>
<pre>
module load git
</pre>
Then you can provide the following line of code to store your GitHub token for future O2 sessions:<br><br>
<pre>
git config --global credential.helper store
</pre>
This process is summarized in the GIF below:<br><br>
<p align="center"><img src="./img/Git_config.gif" width="1000"></p>
Now that we have let Git know to store our GitHub token, we can create one. In order to create a GitHub token, you will need to run:<br><br>
<pre>
usethis::create_github_token()
</pre>
This will take you to a GitHub Webpage. It may prompt you to sign-in to GitHub. From here, you need to name your token and select an expiration date for your token. Then, scroll to the bottom of the page and left-click <kbd>Generate token</kbd>. These step are summarized in the GIF below:<br><br>
<p align="center"><img src="./img/GitHub_token_Part_1.gif" width="1000"></p>
Next, you will want to copy your GitHub Token and go back to RStudio. We need to set our credentials by using the command:<br><br>
<pre>
gitcreds::gitcreds_set()
</pre>
Paste in your copied GitHub token and hit <kbd>Return/Enter</kbd>. It should return:<br><br>
<pre>
-> Adding new credentials...
-> Removing credentials from cache...
-> Done.
</pre>
These steps are summarized in the GIF below:<br><br>
<p align="center"><img src="./img/GitHub_token_Part_2.gif" width="1000"></p>
Now you should have a hidden file called <code>.git-credentials</code> in your Home directory and it should look like:<br><br>
<pre>
https://PersonalAccessToken:YOUR_GITHUB_TOKEN@github.com
</pre>
Where <code>YOUR_GITHUB_TOKEN</code> has been replaced with your GitHub Token.<br><br>
<blockquote>
Note: In order to see hidden file in your file browser on the O2 Portal, you will need to click on the cog in the top right of the <code>Files</code> tab and select <code>Show Hidden Files</code>.
</blockquote>
Next, try pushing your repository to the HBC GitHub again with:<br><br>
<pre>
usethis::use_github(org="hbc",private=TRUE)
</pre>
<hr />
</details>

## Making your first edit

We will now make an edit and push that to GitHub. 

### Editing the `README.md`

We are going to edit the `README.MD` using these steps:

1. Open the `README.md`
2. Change the header from `Guidelines` to the assigned HBC code for the project (`PI_brief_description_hbcXXXX`)
3. Save the changes to `README.md`
4. Select the checkbox in the "Staged" area of the Git tab for the `README.md`
5. Left-click <kbd>Commit</kbd>
6. Add a commit description
7. Left-click commit again
8. Left-click <kbd>Close</kbd> to close the pop-up window for the commit

These steps are summarized in the GIF below:

<p align="center"><img src="./img/Guideline_commit_push.gif" width="1000"></p><br>

### Pushing the edit to GitHub

Now that we have made this commit, we will push it to GitHub using these steps:

1. After completing your commit, left-click <kbd>Push</kbd>
2. Left-click <kbd>Close</kbd> to close the pop-up window for the push
3. Close the additional pop-up window
4. Navigate to the GitHub page for this repository
5. Refresh the page if necessary

You should now see the HBC code as the header to the `README.md` on GitHub. These steps are summarized in the GIF below:

<p align="center"><img src="./img/Guideline_push.gif" width="1000"></p><br>

# Tips for Moving Forward

Now that we've gotten set-up for our project, here are a few last tips to try to make your experience smooth:

- You are welcome to selectively commit and push the parts of these template reports that you would like to have on GitHub.
- Try to avoid editing files directly on GitHub. If you do, it will be important that you `Pull` the repository onto O2 before continuing on with your work on O2. If you forget to do this pull and make commits on O2, you can fix it, but it is beyond the scope of this guide.
- Use the checklist in the `README.md` to help keep track of your progress.


# Note
>These materials have been developed by members of the teaching and platform team at the Harvard Chan Bioinformatics Core (HBC) RRID:SCR_025373. 
