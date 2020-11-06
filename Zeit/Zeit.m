alpha_0=pi/6;
alpha_1=pi/5;

a_0_abs=sin(alpha_0+alpha_1)/sin(alpha_1-2*alpha_0)
a_1_abs=sin(alpha_0+alpha_1)/sin(alpha_0-2*alpha_1)
a_2_abs=1

a_0=a_0_abs*e^(i*alpha_0);
a_1=a_1_abs*e^(i*alpha_1);
a_2=a_2_abs;

% normalization
N=1/sqrt(a_0_abs^2+a_1_abs^2+a_2_abs^2);

% normalized orthogonal row vectors
u_0=[a_0,a_2,a_1]*N;
u_1=[a_1,a_0,a_2]*N;
u_2=[a_2,a_1,a_0]*N;

% check normalization
u_abs_square=u_0*u_0'

% check orthogonality
u_0_u_1=u_0*u_1'
u_1_u_2=u_1*u_2'
u_2_u_3=u_2*u_0'

% unitary matrix
U=[
u_0;
u_1;
u_2;
]

% null matrix
Null=U-U

% unistochastic matrix
P=[
[Null,U.*conj(U)];
[U.*conj(U),Null];
]

% one clock tick
P_square=P*P

% external initial states
psi_0=[1,0,0,0,0,0]';
psi_1=[0,1,0,0,0,0]';
psi_2=[0,0,1,0,0,0]';

% time eigenvalues
t=[0,1,2,0,0,0];

exp_t_0 = t*P_square*psi_0
exp_t_1 = t*P_square*psi_1
exp_t_2 = t*P_square*psi_2

% entropy increase with clock ticks
exp_t_0_2_ticks = t*P_square^2*psi_0
exp_t_1_2_ticks = t*P_square^2*psi_1
exp_t_2_2_ticks = t*P_square^2*psi_2

exp_t_0_3_ticks = t*P_square^3*psi_0
exp_t_1_3_ticks = t*P_square^3*psi_1
exp_t_2_3_ticks = t*P_square^3*psi_2

exp_t_0_4_ticks = t*P_square^4*psi_0
exp_t_1_4_ticks = t*P_square^4*psi_1
exp_t_2_4_ticks = t*P_square^4*psi_2

