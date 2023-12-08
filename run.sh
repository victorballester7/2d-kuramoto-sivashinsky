#add colors to echo's
YELLOW='\033[1;33m'
GREEN='\033[1;32m'
RED='\033[1;31m'
RESET='\033[0m'

echo -e "${YELLOW}Compiling...${RESET}"
# make
cmake -B build -S .
cmake --build build
if [ $? -ne 0 ]; then
    echo -e "${RED}Compilation failed!${RESET}"
    exit 1
fi
echo -e "${GREEN}Compiling done!${RESET}"
echo -e "${YELLOW}Running...${RESET}"
./build/main
if [ $? -ne 0 ]; then
    echo -e "${RED}Running failed!${RESET}"
    exit 1
fi
echo -e "${GREEN}Running done!${RESET}"

# echo -e "${YELLOW}Plotting...${RESET}"
# # python src/plot.py
# julia src/plot.jl
# if [ $? -ne 0 ]; then
#     echo -e "${RED}Plotting failed!${RESET}"
#     exit 1
# fi
# echo -e "${GREEN}Plotting done!${RESET}"