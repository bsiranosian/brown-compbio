# Simulation to answer the Hora and Tempus puzzle!
import math
import random


# can give one a penalty here - number of steps to assemble a watch
horaSteps = 100
tempusSteps = 125
# chance of a call per time step
# vary from 0.0001 to 0.01
calls = [r/10000.0 for r in range(1,101)][::-1]
hora =[]
tempus = []
for c in calls:
	# keep count of watches produced
	horaCount = 0
	tempusCount = 0
	# progress towards current watch 
	horaProgress = 0
	tempusProgress = 0
	# simulate time steps
	for i in range(1000000):
		# is there a call?
		if random.random() < c:
			# if so, reset horaProgress
			horaProgress = 0
			# but only knock 10 parts off Tempus
			if tempusProgress > 10:
				tempusProgress -= 10
			else: tempusProgress = 0

		# no call, each progresses one unit
		horaProgress += 1
		tempusProgress += 1
		
		# check if they've completed a watch
		if horaProgress == horaSteps:
			horaCount +=1
			horaProgress = 0
		if tempusProgress == tempusSteps:
			tempusCount += 1
			tempusProgress = 0
	hora.append(horaCount)
	tempus.append(tempusCount)

plt.plot(calls, hora, label="Hora")
plt.plot(calls, tempus, label="Tempus")
plt.legend()
plt.xlabel('Call Frequency')
plt.ylabel('Watches Produced')
plt.show()