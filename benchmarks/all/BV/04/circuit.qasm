OPENQASM 2.0;
include "qelib1.inc";
qreg qubits[5];

h qubits[0];
h qubits[1];
h qubits[2];
h qubits[3];
h qubits[4];
z qubits[4];
cx qubits[0], qubits[4];
cx qubits[2], qubits[4];
h qubits[0];
h qubits[1];
h qubits[2];
h qubits[3];
h qubits[4];
