#add colors to echo's
YELLOW='\033[1;33m'
BLUE='\033[1;34m'
GREEN='\033[1;32m'
RED='\033[1;31m'
RESET='\033[0m'

# scaled by 100 to avoid floating point errors
NU_MAX=100
NU_MIN=5
NU_STEP=5
# to run the program for several nu's, run it by:
# ./run.sh 1


echo -e "${YELLOW}Compiling...${RESET}"
# make
make
if [ $? -ne 0 ]; then
  echo -e "${RED}Compilation failed!${RESET}"
  exit 1
fi
echo -e "${GREEN}Compilation done!${RESET}"
#
if [ $# -eq 1 ]; then # if there is one argument
  echo -e "${BLUE}Multi run detected.${RESET}"
  read -p "Do you want to remove the old images? (y/n, default: y) " -n 1 -r
  if [[ -z $REPLY ]]; then # if the user just pressed enter
    rm -rf img/energy/*
    echo -e "${YELLOW}Old images removed.${RESET}"
  elif [[ $REPLY =~ ^[Yy]$ ]]; then
    echo 
    rm -rf img/energy/*
    echo -e "${YELLOW}Old images removed.${RESET}"
  else
    echo
  fi
  for ((nu_one=NU_MAX;nu_one>=NU_MIN;nu_one-=NU_STEP)); do
    for ((nu_two=nu_one;nu_two>=NU_MIN;nu_two-=NU_STEP)); do
      echo -e "${YELLOW}Running program for nu1=$nu_one/100 and nu2=$nu_two/100...${RESET}"
      ./bin/main $nu_one $nu_two > /dev/null
      if [ $? -ne 0 ]; then
        echo -e "${RED}Running failed!${RESET}"
        exit 1
      fi
      if [ ! -f "data/tmp_write_E.txt" ]; then
        echo -e "${BLUE}Plotting energy skipped.${RESET}"
        exit 1
      else
        cutoff_time=$(cat data/tmp_write_E.txt)
        rm data/tmp_write_E.txt # remove the tmp file
        python src/plot_energy.py $cutoff_time $nu_one $nu_two > /dev/null
        if [ $? -ne 0 ]; then
          echo -e "${RED}Plotting energy failed!${RESET}"
          continue
        fi
      fi
      echo -e "${GREEN}Program for nu1=$nu_one/100 and nu2=$nu_two/100 done!${RESET}"
    done
  done
  echo -e "${BLUE}Multi run done.${RESET}"
  exit 0
else
  echo -e "${YELLOW}Running...${RESET}"
  ./bin/main
  if [ $? -ne 0 ]; then
    echo -e "${RED}Running failed!${RESET}"
    exit 1
  fi
  echo -e "${GREEN}Running done!${RESET}"
fi

# check if there exists a file called "data/tmp_write_sol.txt"
if [ ! -f "data/tmp_write_sol.txt" ]; then
  echo -e "${BLUE}Animation skipped.${RESET}"
else
  type_anim=$(cat data/tmp_write_sol.txt)
  echo $type_anim
  rm data/tmp_write_sol.txt # remove the tmp file
  echo -e "${YELLOW}Animating...${RESET}"
  python src/animation.py $type_anim
  if [ $? -ne 0 ]; then
    echo -e "${RED}Animating failed!${RESET}"
    exit 1
  fi
  echo -e "${GREEN}Animating done!${RESET}"
fi

# check if there exists a file called "data/tmp_write_E.txt"
if [ ! -f "data/tmp_write_E.txt" ]; then
  echo -e "${BLUE}Plotting energy skipped.${RESET}"
  exit 1
else
  cutoff_time=$(cat data/tmp_write_E.txt)
  rm data/tmp_write_E.txt # remove the tmp file
  echo -e "${YELLOW}Plotting energy...${RESET}"
  python src/plot_energy.py $cutoff_time
  if [ $? -ne 0 ]; then
    echo -e "${RED}Plotting energy failed!${RESET}"
    exit 1
  fi
  echo -e "${GREEN}Plotting energy done!${RESET}"
fi
