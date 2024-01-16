#add colors to echo's
YELLOW='\033[1;33m'
BLUE='\033[1;34m'
GREEN='\033[1;32m'
RED='\033[1;31m'
RESET='\033[0m'

echo -e "${YELLOW}Compiling...${RESET}"
# make
make
if [ $? -ne 0 ]; then
  echo -e "${RED}Compilation failed!${RESET}"
  exit 1
fi
echo -e "${GREEN}Compilation done!${RESET}"
echo -e "${YELLOW}Running...${RESET}"
./bin/main
if [ $? -ne 0 ]; then
  echo -e "${RED}Running failed!${RESET}"
  exit 1
fi
echo -e "${GREEN}Running done!${RESET}"

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
