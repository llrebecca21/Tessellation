arma_spec = function(omega, phi = 0, theta = 0){
  s = rep(NA, length(omega))
  for(i in 1:length(s)){
    z = exp(-2*pi*1i*omega[i])
    s[i] = Mod(1 + crossprod(theta, z^(1:length(theta))))^2 / 
      Mod(1 - crossprod(phi, z^(1:length(phi))))^2
  }
  return(s)
}