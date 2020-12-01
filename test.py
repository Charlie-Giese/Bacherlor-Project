from astropy.visualization import (MinMaxInterval, SqrtStretch,
                                   ImageNormalize)
import numpy as np
import matplotlib.pyplot as plt


data=np.arange(-50, 50, 1).reshape(10,10)

print(data)	

norm = ImageNormalize(data, interval=MinMaxInterval(),
                      stretch=SqrtStretch())

plt.imshow(data, norm=norm)
plt.colorbar()
plt.show()

plt.imshow(data)
plt.colorbar()
plt.show()
