function [xhat_new, P_new, G] = kf(A,B,Bu,C,Q,R,u,y,xhat,P)
	xhat = xhat(:);
	u = u(:);
	y = y(:);

	%�@�\���X�e�b�v
	xhatm = A*xhat + Bu*u;
	Pm = A*P*A' + B*Q*B';

	% �t�B���^�����O�X�e�b�v
	G = Pm*C/(C'*Pm*C+R);

	xhat_new = xhatm + G*(y - C'*xhatm);
	P_new = (eye(size(A)) - G*C')*Pm;
end

