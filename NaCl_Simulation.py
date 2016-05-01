import sys, pygame
pygame.init()

######################################################
################## CONSTANTS USED##################
######################################################

Box_width=200
Box_height=500
Box_depth=300

grey=100,100,100
black=0,0,0
white = 255,255,255
dark=50,50,50

ratio=10  # scaling factor to get size of ions on screen from real size
g=10**-8
G=6.67*(10**-11)

# ionic radii in pm
Actual_Na_Radius=120
Actual_Cl_Radius=170
Na_Radius=Actual_Na_Radius/ratio
Cl_Radius=Actual_Cl_Radius/ratio

Actual_Bond_length=240
Bond_length=Actual_Bond_length/ratio

n_Cl=9
n_Na=7
k=9.00*(10**9)
e=1.60*(10**(-19))


B=e*e*k*(Bond_length**7)/8
mass_Na=23*1.67 * 10**-27
mass_Cl=35.5*1.67 * 10**-27

Frame_Interval=100    # time between frames

######################################################
##################WINDOW FUNCTIONS#################
######################################################

def SetupWindow(height, width, caption):
    screen = pygame.display.set_mode((width, height))
    pygame.display.set_caption(caption)
    return screen

def DrawBackground(screen):
    screen.fill(grey)
    rect=pygame.Rect(Box_width/4, Box_height/4, Box_width/2, Box_height/2)
    pygame.draw.rect(screen,dark,rect)
    pygame.draw.line(screen,black,(0,0),(Box_width/4, Box_height/4))
    pygame.draw.line(screen,black,(3*Box_width/4, 3*Box_height/4),(Box_width, Box_height))
    pygame.draw.line(screen,black,(Box_width/4,3*Box_height/4),(0, Box_height))
    pygame.draw.line(screen,black,(3*Box_width/4,Box_height/4),(Box_width,0))

######################################################
##################ION CLASSES#######################
######################################################

class Ion():
    def __init__(self,_x,_y,_z):
        self.x=_x
        self.y=_y
        self.z=_z
        
        self.fx=0
        self.fy=0
        self.fz=0
        
        self.ax=0
        self.ay=0
        self.az=0
        
        self.vx=0
        self.vy=0
        self.vz=0
        
    def draw(self,screen):
        # get apparent coordinates and radius (take into account perspective)
        c = 1-self.z/(2*Box_depth)
        _x = int((self.x-Box_width/2)*c+Box_width/2)
        _y = int((self.y-Box_height/2)*c+Box_height/2)
        _colour=( self.colour[0] , self.colour[1] , int(255*self.z/Box_depth))
        _r=int(c*self.r)
        pygame.draw.circle(screen, _colour, (int(_x),int(_y)), _r)     # fill
        pygame.draw.circle(screen, black, (int(_x),int(_y)), _r, 1)     # outline
        pygame.draw.circle(screen, white, (int(_x+_r/3),int(_y-_r/3)),int(_r/4))  # glint 

    def ResetForces(self):
        self.fx = 0
        self.fy = 0 #9.81*self.m*g
        self.fz = 0
    
    def bounce(self):
        if self.x>Box_width: #bounce x
            self.x=Box_width
            self.vx=-self.vx     * 0.5
        elif self.x<0:
            self.x=0
            self.vx=-self.vx     * 0.5
        if self.y>Box_height: #bounce y
            self.y=Box_height
            self.vy=-self.vy     * 0.5
        elif self.y<0:
            self.y=0
            self.vy=-self.vy     * 0.5
        if self.z>Box_depth: #bounce z
            self.z=Box_depth
            self.vz=-self.vz    * 0.5
        elif self.z<0:
            self.z=0
            self.vz=-self.vz    * 0.5

    def step(self): # step in time
        self.ax=self.fx/self.m
        self.ay=self.fy/self.m
        self.az=self.fz/self.m
        self.vx+=self.ax*Frame_Interval
        self.vy+=self.ay*Frame_Interval
        self.vz+=self.az*Frame_Interval
        self.x+=self.vx*Frame_Interval
        self.y+=self.vy*Frame_Interval
        self.z+=self.vz*Frame_Interval
        self.bounce()

class Sodium(Ion):
    m = mass_Na
    n = n_Na
    r = Na_Radius
    q = e
    colour = (255,0,0)

class Chlorine(Ion):
    m = mass_Cl
    n = n_Cl
    r = Cl_Radius
    q = -e
    colour = (0,255,0)

def interact(ion1, ion2):
    # calculate the forces acting between the two and add it to their individual net forces
    dx = ion2.x - ion1.x
    dy = ion2.y - ion1.y
    dz = ion2.z - ion1.z
    r2=dx**2+dy**2+dz**2
    if r2==0: # avoid division by 0 error
        dx=1
        r2=1
        
    # calculate forces (negative is attraction)
    Electrostatic = ion1.q*ion2.q*k/r2
    
    r = r2**0.5
    n = (ion1.n + ion2.n) / 2
    repell = 0
    if Electrostatic<0:
        repell = n*B*(r**(-n-1))
    total=(Electrostatic+ repell)/(r)
    ion2.fx+=dx*total
    ion1.fx+=-dx*total
    ion2.fy+=dy*total
    ion1.fy+=-dy*total
    ion2.fz+=dz*total
    ion1.fz+=-dz*total
    
def IonicInteractions(IonList):
    # reset forces on ions:
    for ion in IonList:
        ion.ResetForces()
    # calculate current net forces on ions
    length = len(IonList)
    for index1 in range(length-1):
        for index2 in range(index1+1,length):
            interact( IonList[index1], IonList[index2] )

def SetupIons(width, widthOffset, height, heightOffset, depth, depthOffset):
    IonList = []
    LastIonWasSodium = True
    for _z in range(depth):
        for _y in range(height):
            for _x in range(width):
                x = _x*(Bond_length*1.1)+widthOffset
                y = Box_height-_y*(Bond_length*1.1) - heightOffset
                z = _z*(Bond_length*1.1)+depthOffset
                if(LastIonWasSodium):
                    IonList.append(Chlorine(x,y,z))
                else:
                    IonList.append(Sodium(x,y,z))
                LastIonWasSodium = not LastIonWasSodium
    return IonList

def StepAndDisplay(IonList, screen):
    Cl_list2=[]
    Na_list2=[]
    IonList.sort(key=lambda ion: ion.z, reverse=True)
    for ion in IonList:
        ion.step()
        ion.draw(screen)

def IncrementSpeed(IonList,n):
   for ion in IonList:
      ion.vx *= n
      ion.vy *= n
      ion.vz *= n
def PauseSpeed(IonList):
   for ion in IonList:
      ion.vx = 0
      ion.vy = 0
      ion.vz = 0

if __name__ == "__main__":
    # setup
    screen  = SetupWindow(Box_height, Box_width, "ionic Interactions")
    IonList = SetupIons(5, 20, 5, 20, 4, 20)
    print("Instruction:\n\nhere")
    CarryOn = True
    # main loop
    while CarryOn:
        
        # deal with input
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                pygame.quit()
                CarryOn = False
                break
            elif event.type == pygame.KEYDOWN:
                if event.key == pygame.K_UP:
                    IncrementSpeed(IonList,1.25)
                if event.key == pygame.K_DOWN:
                    IncrementSpeed(IonList,0.8)
                if event.key == pygame.K_RIGHT:
                    PauseSpeed(IonList)

        IonicInteractions(IonList)
        DrawBackground(screen)
        StepAndDisplay(IonList, screen)
        pygame.display.flip()
