Aimer G. Diaz

# Alphafold2 monomer or oligomer prediction

``` bash

module load GCCcore/12.2.0

git config remote.origin.url git@github.com:AimerGDiaz/ProteinStructure
git config --global user.name AimerGDiaz
git config --global user.email "aiagutierrezdi@unal.edu.co"

ssh-keygen -t rsa -C aiagutierrezdi@unal.edu.co

# Generate git.pub in conventional .ssh folder, copy its content into github ssh keys

# Same time write .ssh/config file with 
Host github.com
  HostName github.com
  User git
  IdentityFile ~/.ssh/.git
  IdentitiesOnly yes

#Tesdt key funcion by 
 ssh  -T git@github.com

# Let's start to switching branch, from master to main
git branch main 
git checkout main
git fetch origin
git pull origin main

# Ready to fill 
git add .
git commit -m "Welcome MSU"
git push -u origin main
```
