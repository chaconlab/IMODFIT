################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../libgausscorr.cpp \
../libgaussio.cpp 

OBJS += \
./libgausscorr.o \
./libgaussio.o 

CPP_DEPS += \
./libgausscorr.d \
./libgaussio.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Intel C++ Compiler'
	icpx -O3 -inline-level=2 -I../include -I../../../src -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -c -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


