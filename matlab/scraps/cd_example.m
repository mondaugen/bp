options = optimset('GradObj','on');
x0 = [0;0;0];
x = fminunc(@(x) cd_example_func(x),x0,options)