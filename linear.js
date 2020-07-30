function main() {

	var mathMethods = Object.getOwnPropertyNames(Math);
	for (var i in mathMethods)
		this[mathMethods[i]] = Math[mathMethods[i]];
	
    var ctx = canvas.getContext("2d");
	var w = canvas.width;
	var h = canvas.height;
	
	var ctx_v = canvas_v.getContext("2d");
	var w_v = canvas_v.width;
	var h_v = canvas_v.height;	

	var step = +(document.getElementById('Step').value);

	var particles = [];

	var C0 = 1;
	var m0 = 1;

	var T0 = 0.01;
	
	var R = 1;
	var dt  = step * T0;
	var t = 0;
	var t_bes = 0;
	var tMax = 400 * T0;	
	var fps = 60;	
	var N = +(document.getElementById('Quantity').value);
	var v0 = 0;
	var point = 0;
	var type_print;
	
	var zoom = 8;
	var Chronos = +(document.getElementById('slider').value);
	
	var x = 0;
	var range = 0;
	
	Plus.onclick = function() {
		zoom = zoom + 4;
		if (zoom > 100) {
			zoom = 100
		}
		document.getElementById('Zoom').value = zoom/4;
	}
	
	Minus.onclick = function() {
		zoom = zoom - 4;
		if (zoom < 4) {
			zoom = 4
		}
		document.getElementById('Zoom').value = zoom/4;
	}
	
	New.onclick = function() {
		zoom = 8;
		document.getElementById('Zoom').value = zoom/4;
		N = +(document.getElementById('Quantity').value);
		step = +(document.getElementById('Step').value);
		dt  = step * T0;
		particles = [];
		count();
		point = search_point();
		particles[point].v = +(document.getElementById('v').value);
		v0 = +(document.getElementById('v').value);
		type_print = check('r');		
		ctx_v.clearRect(0,0,w_v,h_v);		
	}
	
	slider.oninput = function(){
		Chronos = +(document.getElementById('slider').value);
		Step.value = slider.value;
	}
	
	r1.onchange = function(){
		check('r');
	}
	
	r2.onchange = function(){
		check('r');
	}

	function search_point(){
		var times;
		if (check('Number') == 0){
			times = 1;
		}
		if (check('Number') == 1) {
			times = Math.round(N/2);
		}
		return times
	}
	
	function check(name){
		var inp = document.getElementsByName(name);
		for (var i = 0; i < inp.length; i++) {
			if (inp[i].type == "radio" && inp[i].checked) {
            type_print = i;
			}
		}
		return type_print
	}	
	
	function control(){
		phys_Bessel();
		draw();
		switch (type_print) {
			case 0: {
				draw_v();
				break
			}
			case 1: {
				draw_v_points();
				break
			}
		}
	}

	function count(){
		for (var i = 1; i < N + 1; i++) {
			var b = [];
			b.x0 = R * i;             
			b.F = 0;   b.v = 0;  b.u = 0; 
			particles[i] = b;               
		}
		
		range = (particles[2].x0 * w_v / N) - (particles[1].x0 * w_v / N);
		
		particles[0] = particles[N]; 
		particles[N+1] = particles[1];
		
		t = 0;
	}
	
	function BesselJ0(x) {
		var ax,z,xx,y,ans,ans1,ans2;
		ax = Math.abs(x)
		if (ax < 8.0) {
			y = x * x;
			ans1 = 57568490574.0 + y * (-13362590354.0 + y * (651619640.7 + y * (-11214424.18 + y * (77392.33017 + y * (-184.9052456)))));
			ans2 = 57568490411.0 + y * (1029532985.0 + y * (9494680.718 + y * (59272.64853 + y * (267.8532712 + y * 1.0))));
			ans = ans1 / ans2;
		} else {
			z = 8.0 / ax;
			y = z * z;
			xx = ax - 0.785398164;
			ans1 = 1.0 + y * (-0.1098628627e-2 + y * (0.2734510407e-4 + y * (-0.2073370639e-5 + y * 0.2093887211e-6)));
			ans2 = -0.1562499995e-1 + y * (0.1430488765e-3 + y * (-0.6911147651e-5 + y * (0.7621095161e-6 - y * 0.934935152e-7)));
			ans = Math.sqrt(0.636619772 / ax) * (Math.cos(xx) * ans1 - z * Math.sin(xx) * ans2);
		}
		return ans;
	}

	function BesselJ1(x) {
		var ax,z,xx,y,ans,ans1,ans2;
		ax = Math.abs(x);
		if (ax < 8.0) {
			y=x*x;
			ans1 = x*(72362614232.0+y*(-7895059235.0+y*(242396853.1+y*(-2972611.439+y*(15704.48260+y*(-30.16036606))))));
			ans2 = 144725228442.0+y*(2300535178.0+y*(18583304.74+y*(99447.43394+y*(376.9991397+y*1.0))));
			ans = ans1/ans2;
		} else {
			z=8.0/ax;
			y=z*z;
			xx=ax-2.356194491;
			ans1=1.0+y*(0.183105e-2+y*(-0.3516396496e-4+y*(0.2457520174e-5+y*(-0.240337019e-6))));
			ans2=0.04687499995+y*(-0.2002690873e-3+y*(0.8449199096e-5+y*(-0.88228987e-6+y*0.105787412e-6)));
			ans=Math.sqrt(0.636619772/ax)*(Math.cos(xx)*ans1-z*Math.sin(xx)*ans2);
			if (x < 0.0) ans = -ans;
		}
		return ans;	
	}
	
	function BesselJn(n,x) {
		var ACC = 40.0;		// Make larger to increase accuracy.
		var BIGNO = 1.0e10;
		var BIGNI = 1.0e-10;
		var j,jsum,m,ax,bj,bjm,bjp,sum,tox,ans;
		ax=Math.abs(x);
		if (ax == 0.0) return 0.0;
		else if (ax > n) {
			tox = 2.0/ax;
			bjm=BesselJ0(ax);
			bj=BesselJ1(ax);
			for (j=1;j<n;j++) {
				bjp=j*tox*bj-bjm;
				bjm=bj;
				bj=bjp;
			}
			ans=bj;
		} else {
			tox=2.0/ax;
			if (Math.sqrt(ACC*n) >= 0)
				m=2*((n + Math.floor(Math.sqrt(ACC*n))) / 2);
			else
				m=2*((n + Math.ceil(Math.sqrt(ACC*n))) / 2);
			jsum=0;
			bjp=ans=sum=0.0;
			bj=1.0;
			for (j=m;j>0;j--) {
				bjm=j*tox*bj-bjp;
				bjp=bj;
				bj=bjm;
				if (Math.abs(bj) > BIGNO) {
					bj *= BIGNI;
					bjp *= BIGNI;
					ans *= BIGNI;
					sum *= BIGNI;
				}
				if (jsum) sum += bj;
				jsum=!jsum;
				if (j == n) ans=bjp;
			}
			sum=2.0*sum-bj;
			ans /= sum;
		}
		return x < 0.0 && (n & 1) ? -ans : ans;
	}

	function Bessel(n,x){
		if (n == 0)
			return BesselJ0(x);
		else if (n == 1)
			return BesselJ1(x);
		else if (n == -1)
			return -BesselJ1(x);
		else if (n > 0)
			return BesselJn(n, x);
		else if (n < 0)
			return Math.pow(-1, -n)*BesselJn(-n, x);
	}

	function phys_Bessel(){
		t_bes = 2*Math.pow(C0/m0, 0.5)*t;
		
		for (var i = 1; i < particles.length-1; i++) {
			particles[i].v_old = particles[i].v;
			particles[i].v = v0 * Bessel(2*(i-point), t_bes);
			particles[i].u += particles[i].v * 2*Math.pow(C0/m0, 0.5)*dt*Chronos/step;
		}

		t += dt*Chronos/step;
		if (t_bes >= 20 * tMax){
			t_bes = 0;
		}
	}

    function draw(){	
		ctx.clearRect(0,0,w,h);

		for (var i=1; i<particles.length - 1; i++){
			ctx.beginPath();
    		ctx.arc((particles[i].x0 + particles[i].u) * w / (N+2), h/2, 500/N, 0, 2*Math.PI);
    		ctx.stroke();
		}
		
	}	
	
	function draw_v(){	
		ctx_v.clearRect(0,0,w_v,h_v);

		ctx_v.beginPath();
		ctx_v.moveTo(0, h_v/2);
		ctx_v.lineTo(w_v, h_v/2);
		ctx_v.strokeStyle = '#000000';
		ctx_v.lineWidth = '1';
		ctx_v.stroke();
		
		
		ctx_v.beginPath();		
		ctx_v.moveTo(particles[1].x0 * w_v / (N+2), h_v/2 - particles[1].v * h_v * zoom);
		for (var i = 2; i<particles.length - 1; i++){	
			ctx_v.lineTo(particles[i].x0 * w_v / (N+2), h_v/2 - particles[i].v * h_v * zoom);
		}
		
		ctx_v.strokeStyle = '#00dd00';
		ctx_v.lineWidth = '2';	
		ctx_v.stroke();
	}
	
	function draw_v_points(){	
		ctx_v.clearRect(0,0,w_v,h_v);

		ctx_v.beginPath();
		ctx_v.moveTo(0, h_v/2);
		ctx_v.lineTo(w_v, h_v/2);
		ctx_v.strokeStyle = '#000000';
		ctx_v.lineWidth = '1';
		ctx_v.stroke();
		
		
		ctx_v.beginPath();	
		ctx_v.fillStyle = 'green';
		ctx_v.moveTo(particles[1].x0 * w_v / (N+2), h_v/2 - particles[1].v_old * h_v * zoom);
		ctx_v.lineTo(particles[1].x0 * w_v / (N+2), h_v/2 - particles[1].v_old * h_v * zoom );
		ctx_v.arc(particles[1].x0 * w_v / (N+2), h_v/2 - particles[1].v * h_v * zoom , 5, 0, 2*Math.PI);
		for (var i = 2; i<particles.length - 1; i++){
			ctx_v.moveTo(particles[i].x0 * w_v / (N+2), h_v/2 - particles[i].v_old * h_v * zoom );
			ctx_v.arc(particles[i].x0 * w_v / (N+2), h_v/2 - particles[i].v * h_v * zoom , 5, 0, 2*Math.PI);
		}
		
		ctx_v.strokeStyle = '#00ff00';
		ctx_v.lineWidth = '0';
		ctx_v.fill();		
		ctx_v.stroke();
	}
	
	count();
	setInterval(control, 1000 / fps);
}