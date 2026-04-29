# Navigate to your folder that you want to push to GitHub
cd /path/to/your/folder

# Initialize a new Git repository in the current folder
# Creates a hidden .git folder that tracks all changes
git init

# Link your local folder to the remote GitHub repository
# 'origin' is just an alias/nickname for the GitHub URL
# This tells Git where to push your code
git remote add origin https://github.com/your-username/your-repository.git

# Stage ALL files in the current folder for commit
# The '.' means everything in the current directory
# Staging = telling Git which files to include in the next commit
git add .

# Create a snapshot/checkpoint of your staged files
# -m flag adds a message describing what this commit contains
# This message appears in your GitHub commit history
git commit -m "Initial commit of scripts"

# Rename the default branch to 'main'
# Older Git versions default to 'master' — GitHub now uses 'main'
# -M forces the rename even if 'main' already exists
git branch -M main

# Store your credentials securely so you only type them once
# Saves username and PAT to ~/.git-credentials file
# Next push will not ask for credentials again
git config --global credential.helper store

# Push your local commits to GitHub
# -u sets 'origin main' as the default upstream
# So next time you can just type 'git push' without specifying origin main
git push -u origin main
# When prompted:
# Username: your-github-username
# Password: paste your new Personal Access Token (NOT your GitHub password)


###############################################################################
# Day 2 — make changes
vim run_batch.py           # edit your file
git add .
git commit -m "Fix GPU loading bug"
# ✅ Local has new commit — GitHub does NOT yet

# Push to sync GitHub
git push
# ✅ Now both identical again

# Day 3 — check status anytime
git status        # shows uncommitted changes
git log --oneline # shows full commit history
git diff          # shows line-by-line changes not yet committed
