################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../pdb2vol.cpp 

OBJS += \
./pdb2vol.o 

CPP_DEPS += \
./pdb2vol.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Intel C++ Compiler'
	icpx -O3 -inline-level=2 -qmkl=sequential -I../../../src -DVOLUME_INCLUDED -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -c -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


