#add colors to echo's
YELLOW='\033[1;33m'
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
echo -e "${GREEN}Compiling done!${RESET}"
echo -e "${YELLOW}Running...${RESET}"
./bin/main
if [ $? -ne 0 ]; then
  echo -e "${RED}Running failed!${RESET}"
  exit 1
fi
echo -e "${GREEN}Running done!${RESET}"

# do not plot if there is an extra argument to ./run.sh
if [ $# -eq 0 ]; then
  echo -e "${YELLOW}Plotting default...${RESET}"
  python src/plot.py 1
  if [ $? -ne 0 ]; then
    echo -e "${RED}Plotting failed!${RESET}"
    exit 1
  fi
else 
  if [ $1 == "1" ]; then
    echo -e "${YELLOW}Plotting surface plot...${RESET}"
    python src/plot.py 1
  elif [ $1 == "2" ]; then
    echo -e "${YELLOW}Plotting contour plot...${RESET}"
    python src/plot.py 2
  else
    echo -e "${RED}Invalid argument! Skipped plotting.${RESET}"
    exit 1
  fi
  echo -e "${GREEN}Plotting done!${RESET}"
fi