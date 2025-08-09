
''' First part: 
      1. A grill in a frequency space is created, Its size is "grid_size". We
      could change the values depending the plot's resolution we wish, but is 
      importat to think in the power of the RAM computer, I'll use 256.
      
      2. The Power spectrum is determinated as k**-3, that's mean that in the
      universe, larger structures are predominant above smaller ones 
      
      
      '''

%reset -f

import numpy as np
import matplotlib.pyplot as plt 

grid_size= 128 # size of the grill 

# spacial frequencies Kx,Ky,Kz in 3 coordenates the function np.fft.fftfreq return normalized values 

kx = np.fft.fftfreq(grid_size)*grid_size
ky = np.fft.fftfreq(grid_size)*grid_size
kz = np.fft.fftfreq(grid_size)*grid_size

# np.meshgrid created a coordenate system with kx,ky,xz ( frequency grill )

kx,ky,kz = np.meshgrid(kx,ky,kz)

# onda number magnitude k 

k = np.sqrt(kx**2 + ky**2 + kz**2)

# Power Spectrum k**-3 (That's mean that larger object in the universe have more influence)

#k[k==0] = 1e-6  # to avoid dividing by zero
#P_k= k**-3
#P_k /= np.max(P_k)

P_k = np.where(k > 0, k**-3, 0)  


#fftshift to move "0,0" to the center

P_k_shifted=np.fft.fftshift(P_k)

# Creating a 3D Gaussian Random field and delta_k 

random_field= np.random.normal(size=(grid_size,grid_size,grid_size))
delta_k=random_field*np.sqrt(P_k) 

delta_real=np.fft.ifftn(delta_k).real

#density field in regular space 

density_field = np.fft.ifftn(delta_k).real 

# ------ Plots -------


plt.figure(figsize=(8,6))
plt.imshow(density_field[:,:,grid_size // 2],cmap='magma',origin='lower')
plt.colorbar(label='density')
plt.title('Perturbations Early Universe')
plt.xlabel('x')
plt.ylabel('y')
plt.show()




''' Plot in Fourier Space : '''

plt.figure(figsize=(8,6))
plt.imshow(np.log10(P_k[:,:,grid_size//2]),cmap='magma',origin='lower')
plt.colorbar(label='Log(P(k))')
plt.title('Power Spectrum')
plt.xlabel('kx')
plt.ylabel('ky')



''' Second part: 
               Solve the poisson equation for gravity 
                (phi_k = Potencial in Fourier)  '''

G = 1 # Gravitational constant (generic value)
rho_bar = 1 # universe density

#phi_k= (4* np.pi * G * rho_bar * delta_k) / k**2 
phi_k = np.where(k > 0, (-4*np.pi* G * rho_bar * delta_k) / k**2+ 1e-6, 0)

phi_real= np.fft.ifftn(phi_k).real

phi_real -= np.mean(phi_real) 

# Plot

plt.figure(figsize=(8,6))
plt.imshow(phi_real[:,:,grid_size//2],cmap='inferno',origin='lower')
plt.colorbar(label='Φ(x)')
plt.title('Potencial gravitacional')
plt.xlabel('x')
plt.ylabel('y')
plt.show()


''' third part : Leapfrog 

      1. We have to compute the aceleration :  a(x) = -∇Φ(x) in Fourier 
      2. From Fourier to real space 
      3. inicialitate velocities 
      4. num_step : It can be change '''


k[k==0] = 1e-6


# Compute aceleration in Fourier 
 
ax_k = np.where(k>0,-1j*kx*phi_k,0) 
ay_k = np.where(k>0,-1j*ky*phi_k,0)
az_k = np.where(k>0,-1j*kz*phi_k,0)


# from Fourier to real space 

ax_real = np.fft.ifftn(ax_k).real
ay_real = np.fft.ifftn(ay_k).real
az_real = np.fft.ifftn(az_k).real 

# velocities star in zero

vx = np.zeros_like(delta_real)
vy = np.zeros_like(delta_real)
vz = np.zeros_like(delta_real)


num_steps = 50
dt= 0.1

for step in range(num_steps):
    # 1. Update the velocities at half-step (v = v + 0.5 a dt)
    vx += 0.5 * ax_real * dt
    vy += 0.5 * ay_real * dt
    vz += 0.5 * az_real * dt
    
    # 2. Update positions (x = x + v dt)
    delta_real += vx * dt
    delta_real += vy * dt
    delta_real += vz * dt
    
    # 3. Recalculate aceleration (as if the density had changed)
    delta_k_new = np.fft.fftn(delta_real)
    phi_k_new = np.where(k > 0, (4 * np.pi * G * rho_bar * delta_k_new) / k**2, 0)
    phi_real_new = np.fft.ifftn(phi_k_new).real  
    
    ax_k_new = np.where(k > 0, -1j * kx * phi_k_new, 0)
    ay_k_new = np.where(k > 0, -1j * ky * phi_k_new, 0)
    az_k_new = np.where(k > 0, -1j * kz * phi_k_new, 0)

    ax_real_new = np.fft.ifftn(ax_k_new).real
    ay_real_new = np.fft.ifftn(ay_k_new).real
    az_real_new = np.fft.ifftn(az_k_new).real

    # 4. Update velocities at full-step (v= v+0.5 a_new dt)
    vx += 0.5 * ax_real_new * dt
    vy += 0.5 * ay_real_new * dt
    vz += 0.5 * az_real_new * dt

# Visualize the final result in a 2D slice
plt.figure(figsize=(8,6))
plt.imshow(delta_real[:, :, grid_size//2], cmap='inferno', origin='lower')
plt.colorbar(label="Density in time")
plt.title("Evolution with gravity (Leapfrog)")
plt.xlabel("x")
plt.ylabel("y")
plt.show()



''' Incorporate the expansion of the universe (Improved)'''

# Compute aceleration in Fourier 
 
ax_k = np.where(k>0,-1j*kx*phi_k,0) 
ay_k = np.where(k>0,-1j*ky*phi_k,0)
az_k = np.where(k>0,-1j*kz*phi_k,0)


# from Fourier to real space 

ax_real = np.fft.ifftn(ax_k).real
ay_real = np.fft.ifftn(ay_k).real
az_real = np.fft.ifftn(az_k).real 

# velocities star in zero

vx = np.zeros_like(delta_real)
vy = np.zeros_like(delta_real)
vz = np.zeros_like(delta_real)


# Parameters
a = 1.0
H0 = 0.05


# Leapfrog with expansion 
for step in range(num_steps):
    a =  (1 + H0 * step * dt)**(2/3) 
    H = H0 * np.sqrt(a)  
    

    # Half-step velocity update with expansion
    vx += 0.5 * ax_real * dt * a  
    vy += 0.5 * ay_real * dt * a  
    vz += 0.5 * az_real * dt * a  

   
    # Cosmological friction 
    vx *= 1- H*dt
    vy *= 1- H*dt
    vz *= 1- H*dt

    # Position update with expansion factor
    delta_real += (vx * dt * a)  # Es la densidad la que se mueve! no en si las particulas
    delta_real += (vy * dt * a) 
    delta_real += (vz * dt * a)  
     

    # Recalculate gravitational potential in Fourier space
    delta_k_new = np.fft.fftn(delta_real)
    phi_k_new = np.where(k > 0, (-4 * np.pi * G * rho_bar * delta_k_new) / (k**2 ), 0)
    print(phi_k_new)
    phi_real_new = np.fft.ifftn(phi_k_new).real  

    # Recalculate acceleration
    ax_k_new = np.where(k > 0, -1j * kx * phi_k_new, 0)
    ay_k_new = np.where(k > 0, -1j * ky * phi_k_new, 0)
    az_k_new = np.where(k > 0, -1j * kz * phi_k_new, 0)

    ax_real_new = np.fft.ifftn(ax_k_new).real
    ay_real_new = np.fft.ifftn(ay_k_new).real
    az_real_new = np.fft.ifftn(az_k_new).real


    # Second half-step with velocities expansion
    vx += 0.5 * ax_real_new * dt * a  
    vy += 0.5 * ay_real_new * dt * a  
    vz += 0.5 * az_real_new * dt * a  

# Visualización en 2D con mejor escala
plt.figure(figsize=(8,6))
plt.imshow(delta_real[:, :, grid_size//2], cmap='inferno', origin='lower')  # Limitar valores extremos
plt.colorbar(label="Final density")
plt.title("Universe with expansion (Fixed Again)")
plt.xlabel("x")
plt.ylabel("y")
plt.show()






