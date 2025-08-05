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
