git status
git log
git remote add origin git@github.com:yourRepository.git
git remote 
git remote -v

git push -u orign main #pushing to main branch
ssh-keygen -t rsa -b 4096 -C "emailAddress@gmail.com" #generates public private key pare in your .ssh/id_rsa.pub
ls ~/.ssh/id_rsa*
cat ~/.ssh/id_rsa*.pub
Put this into github into SSH and GPG keys.
then do:
git push -u orign main #pushing to main branch and all changes pushed to github

