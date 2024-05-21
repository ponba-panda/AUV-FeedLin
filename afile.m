function a = afile(f1,f2,f3,f1d,f2d,f3d,p,phi,q,r,theta,u,v,w)
%AFILE
%    A = AFILE(F1,F2,F3,F1D,F2D,F3D,P,PHI,Q,R,THETA,U,V,W)

%    This function was generated by the Symbolic Math Toolbox version 9.3.
%    12-May-2024 03:17:17

t2 = abs(p);
t3 = abs(q);
t4 = abs(r);
t5 = conj(p);
t6 = cos(phi);
t7 = conj(q);
t8 = conj(r);
t9 = abs(u);
t10 = abs(v);
t11 = cos(theta);
t12 = abs(w);
t13 = conj(u);
t14 = conj(v);
t15 = conj(w);
t16 = sin(phi);
t17 = sin(theta);
t18 = tan(theta);
t19 = f1.*1.000000000000057;
t20 = f1.*-7.139981926858295e-15;
t22 = f1.*1.833133975270958e-15;
t23 = f3.*5.768254396135331e-17;
t24 = f3.*9.296213093854883e-16;
t25 = p.*1.580302784720397e-2;
t26 = q.*1.578973514369978e-2;
t27 = p.*1.060290000000037e+3;
t28 = r.*3.147767789076953e-2;
t29 = q.*-2.007557721793606e-3;
t31 = r.*8.013324904710406e-3;
t32 = f3.*2.171361477651065e-15;
t33 = p.*-8.542600000000675;
t35 = p.*8.586999999999989;
t36 = p.*1.056326674443203e-3;
t37 = q.*-8.586999999999989;
t39 = q.*-8.441790816101286e-3;
t41 = u.*8.079558175371749e-3;
t42 = u.*-6.464602883083601e-2;
t44 = u.*4.068339061975666e-3;
t45 = u.*6.578999490237745e-2;
t46 = p.*4.420190674439551e-3;
t47 = r.*6.878533938992332e-2;
t48 = v.*1.299670777774012e-1;
t49 = w.*6.428530952482388e-2;
t50 = u.*8.542600000000675;
t51 = w.*4.052090694335943e-3;
t52 = w.*6.487042327678548e-2;
t53 = u.*-4.199494918461935e-3;
t55 = u.*4.20486653379904e-3;
t56 = w.*-2.10726614064886;
t58 = p.*1.165149367185861;
t59 = f2.*1.835280624589852e-17;
t60 = r.*5.699675079496197e-1;
t61 = f1.*1.48153079735989e-16;
t62 = q.*1.45698905988319e-1;
t63 = r.*-3.683641833529062e-2;
t65 = p.*-9.673006448961639e-3;
t67 = v.*-1.810979034119594e-2;
t69 = v.*1.159794763754153;
t70 = w.*1.799568424760591e-2;
t71 = v.*3.658516941313077e-2;
t72 = v.*1.840192786276248e-2;
t73 = r.*3.786363251367e-5;
t74 = p.*1.005078337447586e-2;
t75 = w.*1.848353605748874e-2;
t76 = p.*6.465960042325854e-1;
t77 = q.*1.655400000000081e+2;
t78 = u.*-4.812647263489638e-3;
t80 = p.*2.057812373495899e-2;
t81 = p.*4.038914627110604e-5;
t82 = r.*5.052225139634061e-3;
t83 = v.*1.551469605402161e-1;
t84 = f2.*-1.066522753302967e-14;
t86 = v.*1.527758394630965e-4;
t87 = r.*5.168361907175193e-3;
t88 = p.*1.32289283338907e-3;
t89 = v.*4.960570991492592e-3;
t90 = u.*-2.642660000000149e+3;
t92 = u.*1.655400000000081e+2;
t93 = q.*-8.492896240306957e-2;
t95 = r.*1.077399999999943e+1;
t96 = u.*-1.588739457104215e-4;
t98 = r.*8.330161512307703e-5;
t99 = v.*-1.022718570797884e-2;
t102 = w.*1.017414691750599e-2;
t103 = f2.*6.98751912821481e-16;
t104 = r.*2.169986018617404e-2;
t105 = w.*-1.644631437537735e-1;
t107 = q.*5.547485269152687e-3;
t108 = w.*2.074617818167288e-2;
t109 = v.*1.062648251749439e-2;
t110 = w.*1.077499999999964e+1;
t111 = w.*5.522770000000019e+3;
t112 = w.*3.376893993208329e-1;
t113 = u.*-5.445017937384631e-3;
t115 = q.*-5.848799999999756;
t117 = f3.*7.289592919963739e-16;
t118 = q.*1.483799999999974;
t119 = w.*2.730020048496749e-3;
t120 = r.*-1.483799999999974;
t122 = u.*-2.892083448848175;
t124 = q.*-3.864199999999837e+2;
t126 = p.*6.111200000000281;
t127 = u.*1.830972627362553e-1;
t128 = f3.*7.561300199319618e-16;
t129 = f2.*1.891438217392567e-16;
t130 = v.*5.848799999999756;
t132 = r.*3.00547510183367e-3;
t133 = p.*4.807317149915845e-5;
t134 = r.*1.180462193546414e-5;
t135 = v.*1.163221907815916e-2;
t136 = p.*4.094400000000023e+2;
t137 = v.*3.084870000000112e+3;
t138 = q.*-1.636889999999898e+3;
t140 = w.*3.864199999999837e+2;
t141 = q.*2.009226353628769e-1;
t142 = q.*-1.260751668797889e-2;
t144 = w.*6.111200000000281;
t145 = r.*3.16584563657672e-3;
t146 = u.*4.948961243088237e-2;
t147 = v.*-9.815489540701514e-2;
t149 = q.*1.035335019352212e-1;
t150 = u.*6.276075363693023e-3;
t151 = u.*-1.03320000000007e+2;
t153 = v.*4.094400000000023e+2;
t154 = v.*5.129000000000087e+1;
t155 = p.*5.442900000000373e+1;
t156 = p.*1.664544464302284e-3;
t157 = w.*5.129000000000087e+1;
t158 = w.*2.509143019443167e-2;
t159 = r.*1.05696956066895e-1;
t160 = w.*1.03320000000007e+2;
t161 = u.*-1.653254149243628;
t163 = r.*5.442900000000373e+1;
t164 = f2.*2.103263017331757e-16;
t165 = q.*-5.273651243693302e-5;
t167 = p.*3.417434700489158e-3;
t168 = q.*6.792145395852412e-3;
t169 = r.*1.055574639959946e-4;
t170 = f1.*1.069323576631794e-16;
t171 = v.*3.288324231765613e-3;
t172 = p.*-3.47280220025978e-3;
t174 = r.*6.831583685729486e-3;
t175 = v.*3.382170945085102;
t176 = p.*-1.744224225932032e-3;
t178 = f1.*1.083443117011243e-16;
t179 = p.*4.488983132411076e-1;
t180 = u.*-6.877300000000105;
t183 = q.*6.975823768118942e-3;
t185 = q.*1.788596886697117;
t186 = p.*1.127965676993381e-1;
t187 = u.*-1.080874771549432e-1;
t189 = u.*3.381537848214932e-3;
t190 = v.*6.877300000000105;
t191 = r.*9.155499999999884e+2;
t192 = w.*3.455066667689094;
t193 = f2.*9.999999999999364e-1;
t194 = f3.*9.999999999999872e-1;
t195 = q.*3.591297260548953e-3;
t196 = r.*-8.944735076239731e-4;
t198 = p.*1.826824604300971e-3;
t199 = q.*-7.265151547991265e-3;
t201 = v.*2.81562267976818e-2;
t202 = w.*3.495138808153098e-3;
t203 = v.*4.518938472567e-1;
t204 = q.*-1.844492575740158e-3;
t207 = w.*2.257722172700816e-1;
t208 = v.*2.243944459494729e-4;
t209 = q.*-2.417491701184472e-1;
t211 = w.*2.305953886330932e-1;
t212 = r.*7.803299999999581;
t213 = w.*3.671384510515401e-3;
t214 = r.*6.272600000000239e+1;
t216 = r.*2.406331721985364e-4;
t218 = q.*7.778941275231989e-3;
t219 = u.*-7.803299999999581;
t221 = u.*-1.219412271586006e-1;
t223 = w.*-7.507574241134203e-3;
t226 = u.*7.690591039989341e-3;
t227 = v.*6.272600000000239e+1;
t228 = v.*7.744136663763885e-3;
t229 = p.^2;
t230 = q.^2;
t231 = r.^2;
t591 = f1.*3.813830336395979e-4;
t592 = f3.*1.09634965228564e-3;
t593 = f2.*6.255903961139727e-4;
t594 = f3.*5.02896173166604e-5;
t595 = f2.*4.090988023335867e-5;
t596 = f1.*4.090988039000968e-5;
t597 = f2.*4.453513147741838e-5;
t598 = f3.*1.558057199356924e-5;
t599 = f2.*4.570550795193081e-6;
t600 = f1.*3.190541069467161e-6;
t601 = f3.*3.190552021952916e-6;
t605 = f1.*3.40647587366232e-6;
t606 = f2.*1.458290947898744e-6;
t607 = f1.*9.984192481732402e-6;
t608 = f3.*1.980792210023692e-6;
t609 = f1.*1.149055751417578e-6;
t610 = f3.*1.458242124711406e-6;
t611 = f2.*4.013064568116885e-8;
t232 = t3.*3.16822477692391e-2;
t233 = t4.*5.2327e+2;
t235 = t2.*3.379462283679579e-2;
t236 = t2.*2.289151787101079e-3;
t238 = t3.*1.129321031061501e-3;
t239 = t4.*5.736868825515069e-1;
t240 = t11.*1.423280935999937e+2;
t242 = t9.*2.8535841342982e-1;
t243 = t10.*1.140472304954489e-3;
t244 = t12.*1.818131435115952e-2;
t246 = t9.*2.387234833865611e-3;
t247 = t2.*6.72e+2;
t248 = t10.*3.983086975793052e-5;
t249 = t12.*8.109891977169365e-2;
t250 = t17.*1.423280935999918e+2;
t251 = t17.*4.469034445716789;
t252 = t3.*7.7444e+2;
t253 = t9.*7.4822e+2;
t254 = t2.*3.07141013436975e-3;
t255 = t4.*7.630799043069759e-4;
t256 = t17.*3.169320995029377e-1;
t257 = t17.*1.165248371429708e-12;
t258 = t17.*4.058275401740623e-2;
t259 = t17.*7.681942105294646e-11;
t260 = t17.*1.671466028336677e-1;
t261 = t17.*-1.223703439080417e-12;
t263 = t4.*1.669514425420081e-3;
t264 = p.*t5;
t265 = p+t5;
t266 = q.*t6;
t267 = q.*t7;
t268 = q+t7;
t269 = r.*t8;
t270 = r+t8;
t271 = t11.*7.157069849600084e+3;
t272 = t12.*1.82101e+3;
t273 = t13.*u;
t274 = t13+u;
t275 = t3.*-4.84482226366505e-1;
t277 = t12.*3.607042422375244e-3;
t278 = t14.*v;
t279 = t14+v;
t280 = t17.*2.094180000000121e+2;
t281 = t15.*w;
t282 = t15+w;
t283 = t10.*9.9253e+2;
t284 = t9.*3.060959058820362e-2;
t285 = t17.*-2.124722930433663e-1;
t287 = r.*t16;
t288 = t10.*1.546418512077728e-2;
t289 = t17.*-1.689702317442115e-3;
t291 = t17.*-7.157069849599631e+3;
t297 = t6.*t11.*-1.423280935999993e+2;
t299 = t6.*t11.*-8.958494778243932e-3;
t301 = t6.*t11.*-1.811570834848442e-2;
t303 = t6.*t11.*-1.48670353353319e-1;
t306 = t6.*t11.*1.113174484045312e-12;
t307 = t6.*t11.*6.041147510816952e-13;
t308 = t11.*t16.*-2.675248461198126e-4;
t310 = t6.*t17.*1.423280935999937e+2;
t313 = t11.*t16.*1.069486352243648e-12;
t314 = t6.*t11.*1.519398609599934e+3;
t315 = t11.*t16.*1.92920361403157e-2;
t316 = t11.*t16.*-2.277867376537313e-12;
t318 = t11.*t16.*-1.302602978960647;
t320 = t6.*t11.*-1.45660346224803e-12;
t322 = t6.*t11.*2.094179999999897e+2;
t323 = t6.*t17.*1.51939860960003e+3;
t324 = t11.*t16.*-1.519398609600039e+3;
t327 = t6.*t11.*-2.911637480185498e-2;
t329 = t6.*t11.*-9.418452372958379e-1;
t331 = t11.*t16.*-2.094179999999247e+2;
t333 = t6.*t11.*5.958270836301171e-2;
t334 = t16.*t17.*1.51939860960003e+3;
t335 = t11.*t16.*3.167416294863486e-12;
t336 = t11.*t16.*7.15706984960035e+3;
t337 = t11.*t16.*7.027845825233109;
t338 = t11.*t16.*8.913539009973895e-1;
t339 = r.*t6.*t18;
t340 = t11.*t16.*3.048767998535996e-2;
t341 = t16.*t17.*-7.157069849600084e+3;
t343 = q.*t16.*t18;
t407 = t59+t178+t194;
t408 = t128+t129+t170;
t409 = t19+t24+t103;
t410 = t20+t23+t193;
t411 = t22+t84+t117;
t412 = t32+t61+t164;
t419 = t27+t37+t50+t144+t153+t163;
t420 = t35+t92+t120+t130+t138+t140;
t421 = t33+t77+t90+t160+t190+t212;
t422 = t95+t111+t124+t126+t151+t154;
t423 = t115+t136+t137+t157+t180+t214;
t424 = t110+t118+t155+t191+t219+t227;
t602 = -t597;
t603 = -t598;
t604 = -t599;
t612 = -t606;
t613 = -t607;
t614 = -t608;
t615 = -t609;
t616 = -t610;
t617 = -t611;
t822 = t591+t595+t601;
t293 = t247+2.688e+2;
t294 = t233+1.05e+2;
t295 = t253+7.482e+1;
t296 = t252+3.0977e+2;
t305 = t283+6.948e+1;
t311 = t16.*t240;
t312 = t272+7.284e+2;
t326 = t6.*t271;
t344 = -t287;
t345 = 1.0./sqrt(t264);
t346 = 1.0./sqrt(t267);
t347 = 1.0./sqrt(t269);
t348 = 1.0./sqrt(t273);
t349 = 1.0./sqrt(t278);
t350 = 1.0./sqrt(t281);
t413 = t240+t334;
t414 = t271+t323;
t417 = t310+t341;
t418 = p+t339+t343;
t425 = p.*t423.*4.661680209860495e-17;
t426 = p.*t424.*1.835280624589852e-17;
t427 = t423.*u.*-7.289592919963739e-16;
t429 = q.*t424.*2.97996583141159e-16;
t430 = t423.*u.*-7.561300199319618e-16;
t432 = p.*t423.*1.640030243054744e-15;
t433 = t423.*w.*-1.249936544845013e-14;
t435 = p.*t424.*-1.066522753302967e-14;
t437 = p.*t424.*6.98751912821481e-16;
t438 = r.*t423.*1.069323576631794e-16;
t439 = r.*t423.*1.083443117011243e-16;
t440 = t423.*u.*-9.999999999999872e-1;
t442 = p.*t424.*1.891438217392567e-16;
t443 = q.*t424.*1.249936544845013e-14;
t444 = p.*t424.*-2.103263017331757e-16;
t446 = p.*t424.*9.999999999999364e-1;
t447 = q.*t419.*5.768254396135331e-17;
t448 = q.*t419.*9.296213093854883e-16;
t449 = q.*t419.*-2.171361477651065e-15;
t451 = r.*t419.*-1.835280624589852e-17;
t453 = r.*t420.*1.000000000000037;
t454 = p.*t420.*5.768254396135331e-17;
t455 = p.*t420.*9.296213093854883e-16;
t456 = r.*t420.*4.556532371253673e-16;
t457 = r.*t420.*-1.190969582468834e-16;
t459 = r.*t419.*1.066522753302967e-14;
t460 = r.*t420.*-3.980009223755649e-15;
t462 = r.*t419.*-6.98751912821481e-16;
t464 = p.*t420.*-2.171361477651065e-15;
t466 = q.*t419.*7.289592919963739e-16;
t467 = q.*t419.*7.561300199319618e-16;
t468 = r.*t419.*-1.891438217392567e-16;
t470 = r.*t420.*2.97996583141159e-16;
t471 = r.*t419.*2.103263017331757e-16;
t472 = p.*t422.*1.000000000000038;
t473 = q.*t421.*1.000000000000004;
t474 = q.*t422.*-1.000000000000057;
t476 = r.*t421.*1.000000000000038;
t477 = q.*t422.*7.139981926858295e-15;
t478 = p.*t422.*2.904704154411246e-17;
t479 = q.*t422.*-1.833133975270958e-15;
t481 = r.*t421.*2.904704154411246e-17;
t482 = q.*t419.*9.999999999999872e-1;
t483 = r.*t419.*-9.999999999999364e-1;
t485 = t422.*v.*1.000000000000037;
t486 = p.*t420.*7.289592919963739e-16;
t487 = q.*t421.*6.167820657105291e-17;
t488 = t421.*v.*-5.768254396135331e-17;
t490 = p.*t422.*6.215703786013418e-17;
t491 = t422.*v.*4.556532371253673e-16;
t492 = t421.*v.*-9.296213093854883e-16;
t494 = r.*t421.*6.215703786013418e-17;
t495 = p.*t422.*5.103147900608017e-16;
t496 = p.*t422.*4.132340936369327e-15;
t497 = t422.*v.*-1.190969582468834e-16;
t499 = r.*t421.*5.103147900608017e-16;
t500 = p.*t420.*7.561300199319618e-16;
t501 = r.*t421.*4.132340936369327e-15;
t502 = t422.*v.*-3.980009223755649e-15;
t504 = r.*t420.*1.249936544845013e-14;
t505 = t421.*v.*2.171361477651065e-15;
t506 = q.*t422.*1.48153079735989e-16;
t507 = t422.*u.*1.835280624589852e-17;
t508 = t421.*w.*1.835280624589852e-17;
t509 = t422.*v.*2.97996583141159e-16;
t510 = p.*t420.*9.999999999999872e-1;
t511 = q.*t421.*1.649378574341968e-16;
t512 = q.*t421.*1.672374555885555e-16;
t513 = t422.*u.*-1.066522753302967e-14;
t515 = t421.*w.*-1.066522753302967e-14;
t517 = t422.*u.*6.98751912821481e-16;
t518 = t421.*w.*6.98751912821481e-16;
t519 = q.*t421.*4.661680209860495e-17;
t520 = t421.*v.*-7.289592919963739e-16;
t522 = t422.*u.*1.891438217392567e-16;
t523 = t421.*v.*-7.561300199319618e-16;
t525 = t421.*w.*1.891438217392567e-16;
t526 = q.*t421.*1.640030243054744e-15;
t527 = t422.*v.*1.249936544845013e-14;
t528 = p.*t423.*1.000000000000004;
t529 = p.*t422.*-2.114609368560755e-16;
t531 = r.*t421.*-2.114609368560755e-16;
t533 = r.*t423.*1.000000000000057;
t534 = r.*t423.*-7.139981926858295e-15;
t536 = q.*t422.*-1.069323576631794e-16;
t538 = r.*t423.*1.833133975270958e-15;
t539 = q.*t422.*-1.083443117011243e-16;
t541 = t422.*u.*-2.103263017331757e-16;
t543 = t421.*w.*-2.103263017331757e-16;
t545 = t423.*w.*-1.000000000000037;
t547 = p.*t423.*6.167820657105291e-17;
t548 = t423.*u.*-5.768254396135331e-17;
t550 = t423.*u.*-9.296213093854883e-16;
t552 = t423.*w.*-4.556532371253673e-16;
t554 = t422.*u.*9.999999999999364e-1;
t555 = t421.*v.*-9.999999999999872e-1;
t557 = t423.*w.*1.190969582468834e-16;
t558 = t421.*w.*9.999999999999364e-1;
t559 = t423.*w.*3.980009223755649e-15;
t560 = t423.*u.*2.171361477651065e-15;
t561 = r.*t423.*-1.48153079735989e-16;
t563 = q.*t424.*1.000000000000037;
t564 = q.*t424.*4.556532371253673e-16;
t565 = t423.*w.*-2.97996583141159e-16;
t567 = q.*t424.*-1.190969582468834e-16;
t569 = p.*t423.*1.649378574341968e-16;
t570 = p.*t423.*1.672374555885555e-16;
t571 = q.*t424.*-3.980009223755649e-15;
t573 = t45+t56+t62+t67+t88+t196+t232+1.26726536484133e-2;
t574 = t41+t75+t86+t133+t134+t204+t243+7.983639360849332e-5;
t575 = t52+t72+t87+t149+t161+t172+t249+3.243938976815155e-2;
t576 = t31+t65+t89+t105+t122+t141+t288+1.082538142113191e-3;
t577 = t28+t48+t76+t142+t150+t202+t255+1.531205495293681e-4;
t578 = t25+t60+t78+t147+t195+t213+t254+1.2285640537479e-3;
t579 = t70+t93+t98+t109+t146+t167+t239+1.151167134899922e-1;
t580 = t29+t73+t81+t158+t189+t208+t248+2.788277261927612e-6;
t582 = t53+t80+t83+t119+t145+t165+t277+1.442809045781258e-3;
t585 = t69+t104+t108+t113+t186+t199+t263+3.350068122940519e-4;
t586 = t42+t46+t71+t174+t192+t209+t284+3.060877239059895e-3;
t587 = t51+t99+t168+t169+t176+t187+t244+7.272485803693882e-3;
t588 = t63+t96+t107+t198+t201+t223+t236+9.156607148404315e-4;
t590 = t39+t82+t156+t171+t211+t221+t275-1.937891370042253e-1;
t670 = q.*t422.*3.813830336395979e-4;
t671 = r.*t423.*3.813830336395979e-4;
t672 = t422.*v.*9.964608276657059e-4;
t673 = t423.*u.*1.09634965228564e-3;
t675 = p.*t424.*6.255903961139727e-4;
t676 = t423.*w.*9.964608276657059e-4;
t677 = q.*t424.*9.964608276657059e-4;
t678 = q.*t419.*1.09634965228564e-3;
t679 = p.*t422.*3.419504401682718e-4;
t680 = r.*t421.*3.419504401682718e-4;
t681 = p.*t420.*1.09634965228564e-3;
t682 = r.*t419.*6.255903961139727e-4;
t683 = t421.*v.*1.09634965228564e-3;
t684 = r.*t420.*9.964608276657059e-4;
t685 = t422.*u.*6.255903961139727e-4;
t686 = t421.*w.*6.255903961139727e-4;
t688 = t423.*u.*5.02896173166604e-5;
t689 = q.*t421.*1.84402431437929e-4;
t691 = r.*t420.*1.312396551369738e-4;
t693 = r.*t419.*4.090988023335867e-5;
t694 = p.*t422.*1.312396552942353e-4;
t695 = r.*t421.*1.312396552942353e-4;
t696 = p.*t423.*1.84402431437929e-4;
t697 = t422.*v.*1.312396551369738e-4;
t698 = t422.*u.*4.090988023335867e-5;
t699 = t421.*w.*4.090988023335867e-5;
t700 = t423.*w.*1.312396551369738e-4;
t701 = q.*t419.*5.02896173166604e-5;
t702 = q.*t424.*1.312396551369738e-4;
t703 = p.*t420.*5.02896173166604e-5;
t705 = p.*t424.*4.090988023335867e-5;
t706 = t421.*v.*5.02896173166604e-5;
t709 = p.*t423.*9.984191873610243e-6;
t710 = q.*t421.*4.453513165810694e-5;
t711 = t422.*u.*4.453513147741838e-5;
t712 = t421.*w.*4.453513147741838e-5;
t714 = t421.*v.*1.558057199356924e-5;
t715 = t422.*u.*4.570550795193081e-6;
t716 = t421.*w.*4.570550795193081e-6;
t717 = r.*t420.*5.028961773273308e-5;
t718 = p.*t423.*4.453513165810694e-5;
t720 = t423.*u.*1.558057199356924e-5;
t721 = p.*t424.*4.453513147741838e-5;
t722 = t422.*v.*5.028961773273308e-5;
t723 = p.*t422.*1.558056869406205e-5;
t724 = r.*t421.*1.558056869406205e-5;
t725 = p.*t424.*4.570550795193081e-6;
t727 = t423.*w.*5.028961773273308e-5;
t728 = q.*t424.*5.028961773273308e-5;
t729 = q.*t422.*4.090988039000968e-5;
t730 = q.*t421.*9.984191873610243e-6;
t731 = r.*t419.*4.453513147741838e-5;
t732 = q.*t419.*1.558057199356924e-5;
t734 = p.*t420.*1.558057199356924e-5;
t735 = r.*t423.*4.090988039000968e-5;
a = ft_1({f1,f1d,f2,f2d,f3,f3d,p,q,r,t102,t11,t112,t127,t132,t135,t159,t16,t17,t175,t179,t183,t185,t203,t207,t216,t218,t226,t228,t229,t230,t231,t235,t238,t242,t246,t250,t251,t256,t257,t258,t259,t26,t260,t261,t265,t266,t268,t270,t274,t279,t280,t282,t285,t289,t291,t293,t294,t295,t296,t297,t299,t301,t303,t305,t306,t307,t308,t311,t312,t313,t314,t315,t316,t318,t320,t322,t324,t326,t327,t329,t331,t333,t335,t336,t337,t338,t340,t344,t345,t346,t347,t348,t349,t350,t36,t407,t408,t409,t410,t411,t412,t413,t414,t417,t418,t419,t420,t421,t422,t423,t424,t425,t426,t427,t429,t430,t432,t433,t435,t437,t438,t439,t44,t440,t442,t443,t444,t446,t447,t448,t449,t451,t453,t454,t455,t456,t457,t459,t460,t462,t464,t466,t467,t468,t47,t470,t471,t472,t473,t474,t476,t477,t478,t479,t481,t482,t483,t485,t486,t487,t488,t49,t490,t491,t492,t494,t495,t496,t497,t499,t500,t501,t502,t504,t505,t506,t507,t508,t509,t510,t511,t512,t513,t515,t517,t518,t519,t520,t522,t523,t525,t526,t527,t528,t529,t531,t533,t534,t536,t538,t539,t541,t543,t545,t547,t548,t55,t550,t552,t554,t555,t557,t558,t559,t560,t561,t563,t564,t565,t567,t569,t570,t571,t573,t574,t575,t576,t577,t578,t579,t58,t580,t582,t585,t586,t587,t588,t590,t591,t592,t593,t594,t595,t596,t597,t598,t6,t600,t602,t603,t604,t605,t607,t608,t609,t611,t612,t613,t614,t615,t616,t617,t670,t671,t672,t673,t675,t676,t677,t678,t679,t680,t681,t682,t683,t684,t685,t686,t688,t689,t691,t693,t694,t695,t696,t697,t698,t699,t700,t701,t702,t703,t705,t706,t709,t710,t711,t712,t714,t715,t716,t717,t718,t720,t721,t722,t723,t724,t725,t727,t728,t729,t730,t731,t732,t734,t735,t74,t822,u,v,w});
end
function a = ft_1(ct)
[f1,f1d,f2,f2d,f3,f3d,p,q,r,t102,t11,t112,t127,t132,t135,t159,t16,t17,t175,t179,t183,t185,t203,t207,t216,t218,t226,t228,t229,t230,t231,t235,t238,t242,t246,t250,t251,t256,t257,t258,t259,t26,t260,t261,t265,t266,t268,t270,t274,t279,t280,t282,t285,t289,t291,t293,t294,t295,t296,t297,t299,t301,t303,t305,t306,t307,t308,t311,t312,t313,t314,t315,t316,t318,t320,t322,t324,t326,t327,t329,t331,t333,t335,t336,t337,t338,t340,t344,t345,t346,t347,t348,t349,t350,t36,t407,t408,t409,t410,t411,t412,t413,t414,t417,t418,t419,t420,t421,t422,t423,t424,t425,t426,t427,t429,t430,t432,t433,t435,t437,t438,t439,t44,t440,t442,t443,t444,t446,t447,t448,t449,t451,t453,t454,t455,t456,t457,t459,t460,t462,t464,t466,t467,t468,t47,t470,t471,t472,t473,t474,t476,t477,t478,t479,t481,t482,t483,t485,t486,t487,t488,t49,t490,t491,t492,t494,t495,t496,t497,t499,t500,t501,t502,t504,t505,t506,t507,t508,t509,t510,t511,t512,t513,t515,t517,t518,t519,t520,t522,t523,t525,t526,t527,t528,t529,t531,t533,t534,t536,t538,t539,t541,t543,t545,t547,t548,t55,t550,t552,t554,t555,t557,t558,t559,t560,t561,t563,t564,t565,t567,t569,t570,t571,t573,t574,t575,t576,t577,t578,t579,t58,t580,t582,t585,t586,t587,t588,t590,t591,t592,t593,t594,t595,t596,t597,t598,t6,t600,t602,t603,t604,t605,t607,t608,t609,t611,t612,t613,t614,t615,t616,t617,t670,t671,t672,t673,t675,t676,t677,t678,t679,t680,t681,t682,t683,t684,t685,t686,t688,t689,t691,t693,t694,t695,t696,t697,t698,t699,t700,t701,t702,t703,t705,t706,t709,t710,t711,t712,t714,t715,t716,t717,t718,t720,t721,t722,t723,t724,t725,t727,t728,t729,t730,t731,t732,t734,t735,t74,t822,u,v,w] = ct{:};
t736 = r.*t419.*4.570550795193081e-6;
t738 = r.*t420.*4.570548556224273e-6;
t739 = r.*t419.*1.458290947898744e-6;
t740 = r.*t423.*3.190541069467161e-6;
t743 = t423.*u.*3.190552021952916e-6;
t744 = q.*t419.*1.980792210023692e-6;
t745 = t422.*v.*4.570548556224273e-6;
t746 = p.*t420.*1.980792210023692e-6;
t749 = q.*t422.*3.40647587366232e-6;
t750 = q.*t419.*1.458242124711406e-6;
t751 = t422.*u.*1.458290947898744e-6;
t752 = t421.*w.*1.458290947898744e-6;
t753 = q.*t422.*9.984192481732402e-6;
t754 = p.*t422.*1.149055605370427e-6;
t755 = r.*t421.*1.149055605370427e-6;
t756 = t421.*v.*1.980792210023692e-6;
t757 = p.*t420.*1.458242124711406e-6;
t759 = r.*t420.*3.406476377264265e-6;
t760 = t423.*w.*4.570548556224273e-6;
t761 = r.*t420.*4.703663337771495e-7;
t764 = q.*t424.*4.570548556224273e-6;
t765 = r.*t423.*3.40647587366232e-6;
t766 = q.*t422.*1.149055751417578e-6;
t767 = q.*t421.*1.980994375641342e-6;
t769 = t421.*v.*1.458242124711406e-6;
t770 = r.*t423.*9.984192481732402e-6;
t771 = p.*t422.*2.981389479636874e-6;
t772 = r.*t421.*2.981389479636874e-6;
t773 = t422.*v.*3.406476377264265e-6;
t774 = p.*t424.*1.458290947898744e-6;
t775 = t423.*u.*1.980792210023692e-6;
t776 = t422.*v.*4.703663337771495e-7;
t777 = q.*t419.*3.190552021952916e-6;
t779 = q.*t421.*2.981386600636009e-6;
t780 = p.*t420.*3.190552021952916e-6;
t781 = r.*t423.*1.149055751417578e-6;
t782 = p.*t423.*1.980994375641342e-6;
t783 = t423.*u.*1.458242124711406e-6;
t785 = t423.*w.*3.406476377264265e-6;
t786 = t423.*w.*4.703663337771495e-7;
t787 = q.*t422.*3.190541069467161e-6;
t788 = q.*t424.*3.406476377264265e-6;
t790 = p.*t423.*2.981386600636009e-6;
t791 = t421.*v.*3.190552021952916e-6;
t792 = q.*t424.*4.703663337771495e-7;
t795 = p.*t424.*4.013064568116885e-8;
t797 = q.*t421.*4.703756063731425e-7;
t803 = p.*t423.*4.703756063731425e-7;
t807 = r.*t419.*4.013064568116885e-8;
t812 = t422.*u.*4.013064568116885e-8;
t813 = t421.*w.*4.013064568116885e-8;
t818 = p.*t422.*4.013134391157693e-8;
t819 = r.*t421.*4.013134391157693e-8;
t823 = t593+t596+t616;
t824 = t592+t600+t612;
t825 = t594+t604+t605;
t826 = t597+t607+t614;
t827 = t598+t611+t615;
t351 = p.*t293.*1.000000000000037;
t352 = p.*t293.*4.556532371253673e-16;
t353 = p.*t293.*-1.190969582468834e-16;
t355 = p.*t293.*-3.980009223755649e-15;
t357 = p.*t293.*2.97996583141159e-16;
t358 = r.*t294.*-5.768254396135331e-17;
t360 = r.*t294.*-9.296213093854883e-16;
t362 = t295.*u.*-1.000000000000057;
t364 = t295.*u.*7.139981926858295e-15;
t365 = t295.*u.*-1.833133975270958e-15;
t367 = r.*t294.*2.171361477651065e-15;
t368 = q.*t296.*-1.835280624589852e-17;
t370 = t295.*u.*1.48153079735989e-16;
t371 = p.*t293.*1.249936544845013e-14;
t372 = t305.*v.*-1.000000000000038;
t374 = q.*t296.*1.066522753302967e-14;
t375 = t305.*v.*-2.904704154411246e-17;
t377 = q.*t296.*-6.98751912821481e-16;
t379 = t305.*v.*-6.215703786013418e-17;
t381 = t305.*v.*-5.103147900608017e-16;
t383 = r.*t294.*-7.289592919963739e-16;
t385 = t305.*v.*-4.132340936369327e-15;
t387 = t312.*w.*1.000000000000004;
t388 = r.*t294.*-7.561300199319618e-16;
t390 = q.*t296.*-1.891438217392567e-16;
t392 = t312.*w.*6.167820657105291e-17;
t393 = q.*t296.*2.103263017331757e-16;
t394 = t295.*u.*-1.069323576631794e-16;
t396 = r.*t294.*-9.999999999999872e-1;
t398 = t295.*u.*-1.083443117011243e-16;
t400 = q.*t296.*-9.999999999999364e-1;
t402 = t312.*w.*1.649378574341968e-16;
t403 = t312.*w.*1.672374555885555e-16;
t404 = t312.*w.*4.661680209860495e-17;
t405 = t305.*v.*2.114609368560755e-16;
t406 = t312.*w.*1.640030243054744e-15;
t415 = t266+t344;
t416 = t311+t326;
t618 = t295.*u.*3.813830336395979e-4;
t619 = t305.*v.*3.419504401682718e-4;
t620 = r.*t294.*1.09634965228564e-3;
t621 = p.*t293.*9.964608276657059e-4;
t622 = q.*t296.*6.255903961139727e-4;
t624 = t312.*w.*1.84402431437929e-4;
t625 = t305.*v.*1.312396552942353e-4;
t626 = q.*t296.*4.090988023335867e-5;
t629 = r.*t294.*5.02896173166604e-5;
t630 = p.*t293.*1.312396551369738e-4;
t631 = p.*t293.*5.028961773273308e-5;
t634 = t305.*v.*1.558056869406205e-5;
t635 = t295.*u.*4.090988039000968e-5;
t637 = t312.*w.*9.984191873610243e-6;
t638 = q.*t296.*4.453513147741838e-5;
t639 = r.*t294.*1.558057199356924e-5;
t640 = q.*t296.*4.570550795193081e-6;
t641 = t312.*w.*4.453513165810694e-5;
t642 = t295.*u.*9.984192481732402e-6;
t643 = r.*t294.*1.980792210023692e-6;
t644 = p.*t293.*3.406476377264265e-6;
t645 = t305.*v.*1.149055605370427e-6;
t646 = p.*t293.*4.703663337771495e-7;
t647 = r.*t294.*1.458242124711406e-6;
t648 = t295.*u.*1.149055751417578e-6;
t649 = t305.*v.*2.981389479636874e-6;
t650 = t312.*w.*1.980994375641342e-6;
t651 = t312.*w.*2.981386600636009e-6;
t652 = t295.*u.*3.190541069467161e-6;
t653 = r.*t294.*3.190552021952916e-6;
t655 = p.*t293.*4.570548556224273e-6;
t656 = t295.*u.*3.40647587366232e-6;
t657 = q.*t296.*1.458290947898744e-6;
t659 = t312.*w.*4.703756063731425e-7;
t664 = q.*t296.*4.013064568116885e-8;
t668 = t305.*v.*4.013134391157693e-8;
t674 = -t670;
t687 = -t673;
t690 = -t675;
t692 = -t676;
t704 = -t683;
t707 = -t685;
t708 = -t686;
t713 = -t688;
t719 = -t693;
t726 = -t700;
t733 = -t706;
t737 = -t709;
t741 = -t711;
t742 = -t712;
t747 = -t715;
t748 = -t716;
t758 = -t721;
t762 = -t723;
t763 = -t724;
t768 = -t725;
t778 = -t727;
t784 = -t730;
t789 = -t732;
t793 = -t734;
t794 = -t735;
t796 = -t743;
t798 = -t749;
t799 = -t751;
t800 = -t752;
t801 = -t756;
t802 = -t760;
t804 = -t766;
t805 = -t769;
t806 = -t770;
t808 = -t774;
t809 = -t775;
t810 = -t783;
t811 = -t785;
t814 = -t786;
t815 = -t787;
t816 = -t791;
t817 = -t795;
t820 = -t812;
t821 = -t813;
t623 = -t618;
t627 = -t619;
t628 = -t620;
t632 = -t625;
t633 = -t626;
t636 = -t629;
t654 = -t637;
t658 = -t643;
t660 = -t645;
t661 = -t647;
t662 = -t648;
t663 = -t649;
t665 = -t652;
t666 = -t653;
t667 = -t656;
t669 = -t668;
t828 = t250+t307+t324+t355+t368+t381+t396+t398+t404+t425+t426+t439+t440+t451+t460+t482+t495+t499+t502+t507+t508+t510+t519+t539+t555+t559+t571;
t829 = t291+t314+t335+t352+t358+t364+t375+t400+t406+t432+t446+t447+t454+t456+t477+t478+t481+t483+t488+t491+t526+t534+t548+t552+t554+t558+t564;
t830 = t280+t306+t316+t353+t360+t362+t377+t379+t402+t437+t448+t455+t457+t462+t474+t490+t492+t494+t497+t511+t517+t518+t533+t550+t557+t567+t569;
t831 = t259+t313+t322+t357+t365+t374+t383+t387+t405+t427+t429+t435+t459+t466+t470+t473+t479+t486+t509+t513+t515+t520+t528+t529+t531+t538+t565;
t832 = t257+t297+t336+t351+t367+t370+t385+t392+t393+t444+t449+t453+t464+t471+t485+t487+t496+t501+t505+t506+t541+t543+t545+t547+t560+t561+t563;
t833 = t261+t320+t331+t371+t372+t388+t390+t394+t403+t430+t433+t438+t442+t443+t467+t468+t472+t476+t500+t504+t512+t522+t523+t525+t527+t536+t570;
t834 = t258+t303+t337+t621+t632+t636+t640+t659+t667+t672+t677+t684+t692+t694+t695+t701+t703+t713+t733+t736+t747+t748+t765+t768+t797+t798+t803;
t835 = t285+t315+t333+t623+t633+t644+t654+t660+t666+t671+t674+t698+t699+t705+t719+t737+t754+t755+t759+t773+t777+t780+t784+t788+t796+t811+t816;
t836 = t260+t299+t318+t628+t631+t634+t650+t657+t665+t678+t681+t687+t704+t717+t722+t728+t739+t740+t762+t763+t767+t778+t782+t799+t800+t808+t815;
t837 = t289+t301+t338+t627+t630+t639+t651+t662+t664+t679+t680+t691+t697+t702+t714+t720+t726+t779+t781+t789+t790+t793+t804+t807+t817+t820+t821;
t838 = t251+t329+t340+t622+t635+t641+t655+t661+t669+t682+t690+t707+t708+t710+t718+t729+t738+t745+t750+t757+t764+t794+t802+t805+t810+t818+t819;
t839 = t256+t308+t327+t624+t638+t642+t646+t658+t663+t689+t696+t731+t741+t742+t744+t746+t753+t758+t761+t771+t772+t776+t792+t801+t806+t809+t814;
t840 = t822+t835;
t841 = t603+t609+t617+t837;
t842 = t602+t608+t613+t839;
et1 = f1d.*3.813830336395979e-4+f2d.*4.090988023335867e-5+f3d.*3.190552021952916e-6-u.*(t230+t231)-f1.*t573.*4.090988023335867e-5;
et2 = f1.*t574.*(-1.149055605370427e-6)-f2.*t575.*9.984191873610243e-6+f2.*t577.*3.190552021952916e-6+f3.*t576.*1.149055605370427e-6;
et3 = f2.*t578.*3.406476377264265e-6+f2.*t580.*1.149055605370427e-6-f3.*t579.*3.190552021952916e-6+f3.*t582.*9.984191873610243e-6;
et4 = f1.*t585.*(-3.190552021952916e-6)-f1.*t587.*9.984191873610243e-6-f2.*t586.*3.813830336395979e-4-f1.*t588.*3.406476377264265e-6;
et5 = p.*t407.*4.090988023335867e-5-p.*t408.*9.984191873610243e-6-p.*t410.*3.190552021952916e-6-p.*t411.*1.149055605370427e-6;
et6 = p.*t828.*4.090988023335867e-5-p.*t829.*3.190552021952916e-6-p.*t831.*1.149055605370427e-6-p.*t833.*9.984191873610243e-6;
et7 = q.*t407.*3.406476377264265e-6+q.*t409.*9.984191873610243e-6+q.*t411.*3.813830336395979e-4+q.*t412.*3.190552021952916e-6;
et8 = q.*t828.*3.406476377264265e-6+q.*t830.*9.984191873610243e-6+q.*t831.*3.813830336395979e-4-q.*t832.*3.190552021952916e-6-q.*t842.*2.0;
et9 = r.*t408.*3.813830336395979e-4-r.*t409.*1.149055605370427e-6-r.*t410.*3.406476377264265e-6-r.*t412.*4.090988023335867e-5;
et10 = r.*t829.*(-3.406476377264265e-6)-r.*t830.*1.149055605370427e-6+r.*t832.*4.090988023335867e-5+r.*t833.*3.813830336395979e-4-r.*t841.*2.0+t11.*t415.*7.986847213873928e-2;
et11 = t413.*t415.*3.190552021952916e-6-t414.*t415.*4.090988023335867e-5+t415.*t417.*3.406476377264265e-6+t416.*t418.*3.406476377264265e-6;
et12 = t293.*t834.*(-3.406476377264265e-6)-t294.*t836.*3.190552021952916e-6-t295.*t835.*3.813830336395979e-4+t296.*t838.*4.090988023335867e-5;
et13 = t305.*t837.*(-1.149055605370427e-6)+t312.*t839.*9.984191873610243e-6+t590.*t595-t420.*t834.*3.190552021952916e-6-t419.*t836.*4.090988023335867e-5;
et14 = t420.*t836.*3.406476377264265e-6-t422.*t834.*1.149055605370427e-6-t419.*t838.*3.190552021952916e-6+t421.*t836.*1.149055605370427e-6;
et15 = t422.*t835.*4.090988023335867e-5+t423.*t834.*9.984191873610243e-6-t421.*t837.*3.190552021952916e-6-t423.*t835.*3.190552021952916e-6;
et16 = t424.*t834.*(-4.090988023335867e-5)+t421.*t838.*9.984191873610243e-6+t422.*t837.*3.406476377264265e-6+t423.*t836.*3.813830336395979e-4;
et17 = t421.*t839.*(-4.090988023335867e-5)+t422.*t838.*3.813830336395979e-4+t423.*t839.*3.406476377264265e-6-t424.*t838.*3.406476377264265e-6;
et18 = t408.*u.*(-3.190552021952916e-6)-t411.*u.*4.090988023335867e-5-t831.*u.*4.090988023335867e-5-t833.*u.*3.190552021952916e-6;
et19 = t409.*v.*3.190552021952916e-6-t411.*v.*3.406476377264265e-6-t824.*v+t830.*v.*3.190552021952916e-6-t831.*v.*3.406476377264265e-6-t836.*v;
et20 = t408.*w.*(-3.406476377264265e-6)-t409.*w.*4.090988023335867e-5+t823.*w-t830.*w.*4.090988023335867e-5-t833.*w.*3.406476377264265e-6-t838.*w;
et21 = f3.*(t74+t112+t127-t132+t135-t185-t235-1.351784913471832e-2).*3.406476377264265e-6-f3.*(t44+t58-t102+t159+t203-t218-t238-4.517196629718522e-4).*4.090988023335867e-5+f3.*(t47+t49+t175+t179-t183-t226-t246-2.387171022825172e-4).*3.813830336395979e-4+t591.*(t26+t36+t55-t207-t216+t228-t242-2.853507857691472e-2)+p.*q.*v+p.*r.*w;
et22 = t6.*t11.*t418.*-5.088353232777294e-3+t6.*t17.*t415.*2.090869493787761e-3-t11.*t16.*t418.*6.006754565168114e-2+t16.*t17.*t415.*2.406329267654699e-4-p.*t265.*t345.*t825.*1.144576062760793e-3-p.*t265.*t345.*t834.*1.144576062760793e-3-q.*t268.*t346.*t823.*1.584112382396114e-2+q.*t268.*t346.*t838.*1.584112382396114e-2-r.*t270.*t347.*t824.*8.347600782636511e-4-r.*t270.*t347.*t836.*8.347600782636511e-4-t274.*t348.*t822.*u.*1.4267920671491e-1-t274.*t348.*t835.*u.*1.4267920671491e-1+t279.*t349.*t827.*v.*5.702360799991548e-4-t279.*t349.*t837.*v.*5.702360799991548e-4-t282.*t350.*t826.*w.*9.090656621881494e-3+t282.*t350.*t839.*w.*9.090656621881494e-3;
et23 = f1d.*1.149055751417578e-6-f2d.*4.013064568116885e-8-f3d.*1.558057199356924e-5-v.*(t229+t231)+f1.*t573.*4.013064568116885e-8;
et24 = f1.*t574.*(-3.419504401682718e-4)+f2.*t575.*2.981386600636009e-6-f2.*t577.*1.558057199356924e-5+f3.*t576.*3.419504401682718e-4;
et25 = f2.*t578.*1.312396551369738e-4+f2.*t580.*3.419504401682718e-4-f3.*t582.*2.981386600636009e-6+f1.*t585.*1.558057199356924e-5;
et26 = f1.*t587.*2.981386600636009e-6-f2.*t586.*1.149055751417578e-6-f1.*t588.*1.312396551369738e-4-f2.*t590.*4.013064568116885e-8;
et27 = p.*t407.*(-4.013064568116885e-8)+p.*t408.*2.981386600636009e-6+p.*t410.*1.558057199356924e-5-p.*t411.*3.419504401682718e-4;
et28 = p.*t828.*(-4.013064568116885e-8)+p.*t829.*1.558057199356924e-5-p.*t831.*3.419504401682718e-4+p.*t833.*2.981386600636009e-6+p.*t842.*2.0;
et29 = q.*t407.*1.312396551369738e-4-q.*t409.*2.981386600636009e-6+q.*t411.*1.149055751417578e-6-q.*t412.*1.558057199356924e-5;
et30 = q.*t828.*1.312396551369738e-4-q.*t830.*2.981386600636009e-6+q.*t831.*1.149055751417578e-6+q.*t832.*1.558057199356924e-5;
et31 = r.*t408.*1.149055751417578e-6-r.*t409.*3.419504401682718e-4-r.*t410.*1.312396551369738e-4+r.*t412.*4.013064568116885e-8;
et32 = r.*t829.*(-1.312396551369738e-4)-r.*t830.*3.419504401682718e-4-r.*t832.*4.013064568116885e-8+r.*t833.*1.149055751417578e-6+r.*t840.*2.0+t11.*t415.*2.406329573503722e-4;
et33 = t413.*t415.*(-1.558057199356924e-5)+t414.*t415.*4.013064568116885e-8+t415.*t417.*1.312396551369738e-4+t416.*t418.*1.312396551369738e-4;
et34 = t293.*t834.*(-1.312396551369738e-4)+t294.*t836.*1.558057199356924e-5-t295.*t835.*1.149055751417578e-6-t296.*t838.*4.013064568116885e-8;
et35 = t305.*t837.*(-3.419504401682718e-4)-t312.*t839.*2.981386600636009e-6+t579.*t598+t420.*t834.*1.558057199356924e-5+t419.*t836.*4.013064568116885e-8;
et36 = t420.*t836.*1.312396551369738e-4-t422.*t834.*3.419504401682718e-4+t419.*t838.*1.558057199356924e-5+t421.*t836.*3.419504401682718e-4;
et37 = t422.*t835.*(-4.013064568116885e-8)-t423.*t834.*2.981386600636009e-6+t421.*t837.*1.558057199356924e-5+t423.*t835.*1.558057199356924e-5;
et38 = t424.*t834.*4.013064568116885e-8-t421.*t838.*2.981386600636009e-6+t422.*t837.*1.312396551369738e-4+t423.*t836.*1.149055751417578e-6;
et39 = t421.*t839.*4.013064568116885e-8+t422.*t838.*1.149055751417578e-6+t423.*t839.*1.312396551369738e-4-t424.*t838.*1.312396551369738e-4;
et40 = t408.*u.*1.558057199356924e-5+t411.*u.*4.013064568116885e-8+t824.*u+t831.*u.*4.013064568116885e-8+t833.*u.*1.558057199356924e-5+t836.*u;
et41 = t409.*v.*(-1.558057199356924e-5)-t411.*v.*1.312396551369738e-4-t830.*v.*1.558057199356924e-5-t831.*v.*1.312396551369738e-4;
et42 = t408.*w.*(-1.312396551369738e-4)+t409.*w.*4.013064568116885e-8+t825.*w+t830.*w.*4.013064568116885e-8-t833.*w.*1.312396551369738e-4+t834.*w;
et43 = f3.*(t74+t112+t127-t132+t135-t185-t235-1.351784913471832e-2).*1.312396551369738e-4+f3.*(t44+t58-t102+t159+t203-t218-t238-4.517196629718522e-4).*4.013064568116885e-8+f3.*(t47+t49+t175+t179-t183-t226-t246-2.387171022825172e-4).*1.149055751417578e-6+t609.*(t26+t36+t55-t207-t216+t228-t242-2.853507857691472e-2)+p.*q.*u;
et44 = q.*r.*w-t6.*t11.*t418.*4.793747785535862e-2-t6.*t17.*t415.*6.243560191320069e-4-t11.*t16.*t418.*5.633815718816876e-4+t16.*t17.*t415.*7.161057727916089e-2-p.*t265.*t345.*t825.*4.40965241260232e-2-p.*t265.*t345.*t834.*4.40965241260232e-2+q.*t268.*t346.*t823.*1.55393886206622e-5-q.*t268.*t346.*t838.*1.55393886206622e-5+r.*t270.*t347.*t824.*4.076422953537489e-3+r.*t270.*t347.*t836.*4.076422953537489e-3-t274.*t348.*t822.*u.*4.298732471628301e-4-t274.*t348.*t835.*u.*4.298732471628301e-4+t279.*t349.*t827.*v.*1.696980351901074e-1-t279.*t349.*t837.*v.*1.696980351901074e-1+t282.*t350.*t826.*w.*2.714567406812089e-3-t282.*t350.*t839.*w.*2.714567406812089e-3;
et45 = f1d.*9.984192481732402e-6+f2d.*4.453513147741838e-5-f3d.*1.980792210023692e-6-w.*(t229+t230)-f1.*t573.*4.453513147741838e-5;
et46 = f1.*t574.*2.981389479636874e-6-f2.*t575.*1.84402431437929e-4-f2.*t577.*1.980792210023692e-6-f3.*t576.*2.981389479636874e-6;
et47 = f2.*t578.*(-4.703663337771495e-7)-f2.*t580.*2.981389479636874e-6+f3.*t582.*1.84402431437929e-4+f1.*t585.*1.980792210023692e-6;
et48 = f1.*t587.*(-1.84402431437929e-4)-f2.*t586.*9.984192481732402e-6+f1.*t588.*4.703663337771495e-7+p.*t407.*4.453513147741838e-5;
et49 = p.*t408.*(-1.84402431437929e-4)+p.*t410.*1.980792210023692e-6+p.*t411.*2.981389479636874e-6+p.*t828.*4.453513147741838e-5;
et50 = p.*t829.*1.980792210023692e-6+p.*t831.*2.981389479636874e-6-p.*t833.*1.84402431437929e-4+p.*t841.*2.0-q.*t407.*4.703663337771495e-7;
et51 = q.*t409.*1.84402431437929e-4+q.*t411.*9.984192481732402e-6-q.*t412.*1.980792210023692e-6-q.*t828.*4.703663337771495e-7;
et52 = q.*t830.*1.84402431437929e-4+q.*t831.*9.984192481732402e-6+q.*t832.*1.980792210023692e-6-q.*t840.*2.0+r.*t408.*9.984192481732402e-6;
et53 = r.*t409.*2.981389479636874e-6+r.*t410.*4.703663337771495e-7-r.*t412.*4.453513147741838e-5+r.*t829.*4.703663337771495e-7;
et54 = r.*t830.*2.981389479636874e-6+r.*t832.*4.453513147741838e-5+r.*t833.*9.984192481732402e-6+t11.*t415.*2.090869621139487e-3-t413.*t415.*1.980792210023692e-6;
et55 = t414.*t415.*(-4.453513147741838e-5)-t415.*t417.*4.703663337771495e-7-t416.*t418.*4.703663337771495e-7+t293.*t834.*4.703663337771495e-7;
et56 = t294.*t836.*1.980792210023692e-6-t295.*t835.*9.984192481732402e-6+t296.*t838.*4.453513147741838e-5+t305.*t837.*2.981389479636874e-6;
et57 = t312.*t839.*1.84402431437929e-4+t579.*t608+t590.*t597+t420.*t834.*1.980792210023692e-6-t419.*t836.*4.453513147741838e-5-t420.*t836.*4.703663337771495e-7;
et58 = t422.*t834.*2.981389479636874e-6+t419.*t838.*1.980792210023692e-6-t421.*t836.*2.981389479636874e-6+t422.*t835.*4.453513147741838e-5;
et59 = t423.*t834.*1.84402431437929e-4+t421.*t837.*1.980792210023692e-6+t423.*t835.*1.980792210023692e-6-t424.*t834.*4.453513147741838e-5;
et60 = t421.*t838.*1.84402431437929e-4-t422.*t837.*4.703663337771495e-7+t423.*t836.*9.984192481732402e-6-t421.*t839.*4.453513147741838e-5;
et61 = t422.*t838.*9.984192481732402e-6-t423.*t839.*4.703663337771495e-7+t424.*t838.*4.703663337771495e-7+t408.*u.*1.980792210023692e-6;
et62 = t411.*u.*(-4.453513147741838e-5)-t823.*u-t831.*u.*4.453513147741838e-5+t833.*u.*1.980792210023692e-6+t838.*u-t409.*v.*1.980792210023692e-6;
et63 = t411.*v.*4.703663337771495e-7-t825.*v-t830.*v.*1.980792210023692e-6+t831.*v.*4.703663337771495e-7-t834.*v+t408.*w.*4.703663337771495e-7;
et64 = t409.*w.*(-4.453513147741838e-5)-t830.*w.*4.453513147741838e-5+t833.*w.*4.703663337771495e-7-f3.*(t74+t112+t127-t132+t135-t185-t235-1.351784913471832e-2).*4.703663337771495e-7;
et65 = f3.*(t44+t58-t102+t159+t203-t218-t238-4.517196629718522e-4).*(-4.453513147741838e-5)+f3.*(t47+t49+t175+t179-t183-t226-t246-2.387171022825172e-4).*9.984192481732402e-6+t607.*(t26+t36+t55-t207-t216+t228-t242-2.853507857691472e-2)+p.*r.*u+q.*r.*v+t6.*t11.*t418.*3.633969551863179e-3+t6.*t17.*t415.*3.861718838686916e-2-t11.*t16.*t418.*2.904942845827486e-2-t16.*t17.*t415.*6.243566220466101e-4+p.*t265.*t345.*t825.*1.580430881491222e-4+p.*t265.*t345.*t834.*1.580430881491222e-4-q.*t268.*t346.*t823.*1.724489361068595e-2+q.*t268.*t346.*t838.*1.724489361068595e-2+r.*t270.*t347.*t824.*5.182445698695488e-4;
et66 = r.*t270.*t347.*t836.*5.182445698695488e-4-t274.*t348.*t822.*u.*3.735186249340909e-3-t274.*t348.*t835.*u.*3.735186249340909e-3-t279.*t349.*t827.*v.*1.479559250111993e-3+t279.*t349.*t837.*v.*1.479559250111993e-3-t282.*t350.*t826.*w.*1.678993358363915e-1+t282.*t350.*t839.*w.*1.678993358363915e-1;
a = [et1+et2+et3+et4+et5+et6+et7+et8+et9+et10+et11+et12+et13+et14+et15+et16+et17+et18+et19+et20+et21+et22;et23+et24+et25+et26+et27+et28+et29+et30+et31+et32+et33+et34+et35+et36+et37+et38+et39+et40+et41+et42+et43+et44;et45+et46+et47+et48+et49+et50+et51+et52+et53+et54+et55+et56+et57+et58+et59+et60+et61+et62+et63+et64+et65+et66];
end
