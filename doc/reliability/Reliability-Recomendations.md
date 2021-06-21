# Reliability Recomendations

The following details recomendations and analysis of the reliability of the initial architecture of the Lunar Gateway amateur radio payload.

An initial risk analysis of the architecture was performed as the basis of these recomendations, the risk analysis can be found here: [Risk Analysis](https://github.com/phase4space/payload-dmt/blob/master/doc/reliability/Risk-Analysis.md)

The system architecture is shown below with hot (red) and cold (blue) redundant systems.

![System Architecture](https://github.com/phase4space/payload-dmt/blob/master/doc/reliability/diagrams/Redundancy-Zones-ARExV3.png "System Architecture")

The feasibility of adding these redundant systems is dependant on the available volume and the architecture will evolve over time.  The following items serve as a starting point for the discussion.

## Antenna and Amplifiers

The antennas, cabling and amplifiers are a single point of failure for the primary mission objective. Cabling and connectors are common failure points due to the high mechanical stresses during launch and temperature cycling. The antennas are exposed to the space environment and power amplifiers are subjected to stress from high power, voltage and temperature.  Therefore these systems are of primary concern.

To mitigate the risk a power splitting and combining topology is recomended.  The C/X band patch antennas will almost certainly be an array with multiple elements. Therefore it is proposed that multiple parallel paths are formed of array elements, cabling, connectors and amplifiers. One side of the power combination/splitting is performed in free-space and the other end is performed in passive microwave power splitter/combiner topologies which are very high reliability.

![Antenna Structure](https://github.com/phase4space/payload-dmt/blob/master/doc/reliability/diagrams/antenna_power_combine.svg "Antenna Structure")

This provides passive redundancy - if a failure occurs the other parallel paths are already active. The topology has the added benefit that the power amplifiers in the transmit paths can be lower powered. In the receive path the signal to noise is improved as the signal in each LNA is correlated but the noise in uncorrelated.

## Cold Redundant Digital Processing

The digital processing system is complex and as such is at risk of failure which would lead to loss of the primary mission objective. A cold redundant processing board can be used that will be automatically switched over if a problem occurs.

By placing the ADC/DACs on the RF board the digital signals can be easily switched or combined to the primary/secondary digital processing board.

## TT&C Control

Using a RadHard TT&C processor as a system controller provides reliable control of the system. A hardware RadHard watchdog can be used to automatically switch between two TT&C processors if one fails.  Ideally two TT&C buses will be used to protect against a failed transceiver bringing the bus (and whole system) down.  The TT&C controllers can be used to check the state of the main digital processing and switch between primary and secondary boards.

![State Transfer](https://github.com/phase4space/payload-dmt/blob/master/doc/reliability/diagrams/state_transfer.svg "State Transfer")


## Clocking and Frequency Synthesis

Clock generation has well known reliability issues in space missions.  Therefore, a simple cold redundant clocking circuit can be used as shown below.

![Redundnacy](https://github.com/phase4space/payload-dmt/blob/master/doc/reliability/diagrams/oscillator_cold_spare.svg "Clock Redundnacy")

The circuit will monitor the output of the primary oscillator - if it stops operating it will be disabled and the secondary oscillator automatically enabled to maintain system operation.

The synthesis of the Rx and Tx RF section could also use a common circuit with a frequency double or divider.  The advantage being that the area that would be used for the secondary LO synthesiser can be used for a cold redundant synthesiser.

## Power

Power is the bedrock of a reliable system.  Therefore, each power bus should be fitted with an LCL (Latching Current Limiter) which detects and protects against latch-up in downstream circuits. Ideally, parallel buses and no-single point bus shorting failure design should be used to protect against power bus failure.

The harnessing should use parallel pins and connectors wherever possible.
