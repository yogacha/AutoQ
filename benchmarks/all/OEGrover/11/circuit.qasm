OPENQASM 2.0;
include "qelib1.inc";
qreg qubits[22];

x qubits[0];
x qubits[2];
x qubits[4];
x qubits[6];
x qubits[8];
x qubits[10];
ccx qubits[0], qubits[1], qubits[12];
ccx qubits[2], qubits[12], qubits[13];
ccx qubits[3], qubits[13], qubits[14];
ccx qubits[4], qubits[14], qubits[15];
ccx qubits[5], qubits[15], qubits[16];
ccx qubits[6], qubits[16], qubits[17];
ccx qubits[7], qubits[17], qubits[18];
ccx qubits[8], qubits[18], qubits[19];
ccx qubits[9], qubits[19], qubits[20];
cz qubits[20], qubits[10];
ccx qubits[9], qubits[19], qubits[20];
ccx qubits[8], qubits[18], qubits[19];
ccx qubits[7], qubits[17], qubits[18];
ccx qubits[6], qubits[16], qubits[17];
ccx qubits[5], qubits[15], qubits[16];
ccx qubits[4], qubits[14], qubits[15];
ccx qubits[3], qubits[13], qubits[14];
ccx qubits[2], qubits[12], qubits[13];
ccx qubits[0], qubits[1], qubits[12];
x qubits[0];
x qubits[2];
x qubits[4];
x qubits[6];
x qubits[8];
x qubits[10];
h qubits[0];
h qubits[1];
h qubits[2];
h qubits[3];
h qubits[4];
h qubits[5];
h qubits[6];
h qubits[7];
h qubits[8];
h qubits[9];
h qubits[10];
x qubits[0];
x qubits[1];
x qubits[2];
x qubits[3];
x qubits[4];
x qubits[5];
x qubits[6];
x qubits[7];
x qubits[8];
x qubits[9];
x qubits[10];
ccx qubits[0], qubits[1], qubits[12];
ccx qubits[2], qubits[12], qubits[13];
ccx qubits[3], qubits[13], qubits[14];
ccx qubits[4], qubits[14], qubits[15];
ccx qubits[5], qubits[15], qubits[16];
ccx qubits[6], qubits[16], qubits[17];
ccx qubits[7], qubits[17], qubits[18];
ccx qubits[8], qubits[18], qubits[19];
ccx qubits[9], qubits[19], qubits[20];
cz qubits[20], qubits[10];
ccx qubits[9], qubits[19], qubits[20];
ccx qubits[8], qubits[18], qubits[19];
ccx qubits[7], qubits[17], qubits[18];
ccx qubits[6], qubits[16], qubits[17];
ccx qubits[5], qubits[15], qubits[16];
ccx qubits[4], qubits[14], qubits[15];
ccx qubits[3], qubits[13], qubits[14];
ccx qubits[2], qubits[12], qubits[13];
ccx qubits[0], qubits[1], qubits[12];
x qubits[0];
x qubits[1];
x qubits[2];
x qubits[3];
x qubits[4];
x qubits[5];
x qubits[6];
x qubits[7];
x qubits[8];
x qubits[9];
x qubits[10];
h qubits[0];
h qubits[1];
h qubits[2];
h qubits[3];
h qubits[4];
h qubits[5];
h qubits[6];
h qubits[7];
h qubits[8];
h qubits[9];
h qubits[10];
z qubits[11];
