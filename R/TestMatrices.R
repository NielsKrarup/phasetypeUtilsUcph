
# Test Matrices -----------------------------------------------------------
listS <- list(

  S1 = matrix(c(-10, 2, 3,
                1, -12, 10,
                6, 1,  -8) , byrow = T, ncol = 3),


  S2 = matrix(c(-4, 1, 1, 1,
                1, -5, 1, 1,
                1, 1, -6, 1,
                1, 1, 1, -7), byrow = T, ncol = 4),


  S3 = matrix(c(-4, 3,
                1, -5), byrow = T, ncol = 2),

  S4 = matrix(c(-7, 1, 1, 1,2,
                1, -7, 1, 1,2,
                1, 1, -10, 1,1,
                1, 1, 1, -8,3,
                3,3,2,1,-10), byrow = T, ncol = 5),

  S5 = matrix(c(-5,1,1,
                3,-4,1,
                3,2,-7), byrow = T, ncol = 3),


  S6 = matrix(c(-4,2,1,
                3,-6,1,
                3,1,-7), byrow = T, ncol = 3))

listAlpha <- list(
  pi1 = c(1,4,5)/10,
  pi2 = c(1,1,1,1)/4,
  pi3 = c(0.1,0.9),
  pi4 = c(1,1,1,1,1)/5,
  pi5 = c(1,2,3)/6,
  pi6 = c(2,3,5)/10)


# listAlpha = list(pi1,pi2,pi3,pi4,pi5,pi6)
# listS = list(S1,S2,S3,S4,S5,S6)

# listStmp = list(S3,S1,S2,S4)
# listAtmp = list(pi3,pi1,pi2,pi4)
