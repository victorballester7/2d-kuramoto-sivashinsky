#add colors to echo's
YELLOW='\033[1;33m'
BLUE='\033[1;34m'
GREEN='\033[1;32m'
RED='\033[1;31m'
RESET='\033[0m'

# scaled by 100 to avoid floating point errors
# set of variables for complete multi run
NU_MAX=100
NU_MIN=5
NU_STEP=5
# set of variables for multi run with nu1=nu2
NU_MAX_SAME=10
NU_MIN_SAME=1
NU_STEP_SAME=1

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
  # check if there is a folder called "img/energy"
  if [ ! -d "img/energy" ]; then
    mkdir img/energy
  else
    read -p "Old images detected. Do you want to remove them? (y/n, default: y) " -n 1 -r
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
  fi
  read -p "Do you want to run the program for nu1 = nu2? (y/n, default: y) " -n 1 -r
  if [[ -z $REPLY ]]; then # if the user just pressed enter
    SAME=1
  elif [[ $REPLY =~ ^[Yy]$ ]]; then
    echo
    SAME=1
  else
    echo
    SAME=0
  fi
  if [ $SAME -eq 1 ]; then
    for ((nu=NU_MAX_SAME;nu>=NU_MIN_SAME;nu-=NU_STEP_SAME)); do
      echo -e "${YELLOW}Running program for nu1=$nu/100 and nu2=$nu/100...${RESET}"
      ./bin/main $nu $nu > /dev/null
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
        python src/plot_energy.py $cutoff_time $nu $nu > /dev/null
        if [ $? -ne 0 ]; then
          echo -e "${RED}Plotting energy failed!${RESET}"
          continue
        fi
      fi
      echo -e "${GREEN}Program for nu1=$nu/100 and nu2=$nu/100 done!${RESET}"
    done
  else
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
  fi
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
  rm data/tmp_write_sol.txt # remove the tmp file
  cutoff_time=$(cat data/tmp_write_E.txt)
  echo -e "${YELLOW}Animating...${RESET}"
  python src/animation.py $type_anim $cutoff_time
  if [ $? -ne 0 ]; then
    echo -e "${RED}Animating failed!${RESET}"
    exit 1
  fi
  echo -e "${GREEN}Animating done!${RESET}"
fi

# check if there exists a file called "data/tmp_write_E.txt"
if [ ! -f "data/tmp_write_E.txt" ]; then
  echo -e "${BLUE}Plotting energy skipped.${RESET}"
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

if [ ! -f "data/tmp_write_slice_sol.txt" ]; then
  echo -e "${BLUE}Plotting slice skipped.${RESET}"
else
  args=$(cat data/tmp_write_slice_sol.txt)
  rm data/tmp_write_slice_sol.txt # remove the tmp file
  echo -e "${YELLOW}Plotting slice...${RESET}"
  python src/plot_solution.py $args
  if [ $? -ne 0 ]; then
    echo -e "${RED}Plotting slice failed!${RESET}"
    exit 1
  fi
  echo -e "${GREEN}Plotting slice done!${RESET}"
fi

if [ ! -f "data/tmp_write_slice_freq.txt" ]; then
  echo -e "${BLUE}Plotting slice frequency skipped.${RESET}"
else
  args=$(cat data/tmp_write_slice_freq.txt)
  rm data/tmp_write_slice_freq.txt # remove the tmp file
  echo -e "${YELLOW}Plotting slice frequency...${RESET}"
  python src/plot_freq.py $args
  if [ $? -ne 0 ]; then
    echo -e "${RED}Plotting slice frequency failed!${RESET}"
    exit 1
  fi
  echo -e "${GREEN}Plotting slice frequency done!${RESET}"
fi

if [ ! -f "data/tmp_estimate_period.txt" ]; then
  echo -e "${BLUE}Estimating period skipped.${RESET}"
else
  rm data/tmp_estimate_period.txt # remove the tmp file
  echo -e "${YELLOW}Estimating period...${RESET}"
  python src/find_period.py $cutoff_time
  if [ $? -ne 0 ]; then
    echo -e "${RED}Estimating period failed!${RESET}"
    exit 1
  fi
  echo -e "${GREEN}Estimating period done!${RESET}"
fi

if [ ! -f "data/tmp_estimate_period_burst.txt" ]; then
  echo -e "${BLUE}Estimating period for burst skipped.${RESET}"
else
  rm data/tmp_estimate_period_burst.txt # remove the tmp file
  echo -e "${YELLOW}Estimating period for burst...${RESET}"
  python src/find_period_burst.py $cutoff_time
  if [ $? -ne 0 ]; then
    echo -e "${RED}Estimating period for burst failed!${RESET}"
    exit 1
  fi
  echo -e "${GREEN}Estimating period for burst done!${RESET}"
fi