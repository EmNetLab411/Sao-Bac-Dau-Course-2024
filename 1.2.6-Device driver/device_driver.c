# include <linux/module.h>
# include <linux/init.h>
# include <linux/slab.h>
# include <linux/i2c.h>
# include <linux/delay.h>
# include <linux/kernel.h>

#define I2C_BUS_AVAILABLE   (          1 ) // Bus I2C khả dụng trên Raspberry Pi
#define SLAVE_DEVICE_NAME   ( "PCA9685" ) // Tên thiết bị 
#define PCA9685_SLAVE_ADDR  (0x40) // Địa chỉ I2C Bus của PCA9685

/*
Định nghĩa các thanh ghi
Sinh viên tra cứu trong tài liệu PCA9685.pdf mục 7.3 (Register definitions) để hoàn thiện phần này
*/

#define SUBADR1          0x02
#define SUBADR2          0x03
#define SUBADR3          0x04
#define MODE1            0x00
#define MODE2            0x01
#define PRESCALE         
#define LED0_ON_L        
#define LED0_ON_H        
#define LED0_OFF_L       
#define LED0_OFF_H       
#define ALLLED_ON_L      
#define ALLLED_ON_H      
#define ALLLED_OFF_L     
#define ALLLED_OFF_H     

#define PWM_I2C_Addr 	 0x40 // Địa chỉ I2C Bus của PCA9685
#define PWM_I2C_Hz 	 50

static struct i2c_adapter *PCA9685_i2c_adapter = NULL;   // I2C Adapter Structure
static struct i2c_client  *PCA9685_i2c_client = NULL;   // I2C Cient Structure (PCA9685)

/*
Những hàm quan trọng để điều khiển cánh tay robot
*/
static void PCA9685_setPWM(uint8_t channel, uint16_t on, uint16_t off);
static void PCA9685_setServoPulse(uint8_t channel, uint16_t value);
static void PCA9685_Set_Rotation_Angle(uint8_t channel, uint8_t Angle);

/*
** Hàm viết dữ liệu vào thanh ghi trên PCA9685 qua I2C. 
** 
** Argument:
**	  reg  -> Thanh ghi truy cập tới
**	  value-> Giá trị 
**    buff -> Byte ghi vào slave 
**    len  -> Số 
*/
static int I2C_Write(unsigned char reg, unsigned char value)
{

	unsigned char buf[2];
	buf[0] = reg;
	buf[1] = value;
    /*
    ** Sending Start condition, Slave address with R/W bit, 
    ** ACK/NACK and Stop condtions will be handled internally.
    */
	int ret = i2c_master_send(PCA9685_i2c_client, buf, 2);

	return ret; // Tra ve -1 neu loi,neu khong thi tra ve so byte duoc viet.
}

/*
** Ham nay doc 1 byte data tu I2C client
**
**  Arguments:
**      out_buff -> buffer wherer the data to be copied
**      len      -> Length of the data to be read
** 
*/
static int I2C_Read(unsigned char reg)
{
    /*
    ** Sending Start condition, Slave address with R/W bit, 
    ** ACK/NACK and Stop condtions will be handled internally.
    */
	int ret = i2c_smbus_read_byte_data(PCA9685_i2c_client, reg);

	return ret; // Tra ve -1 neu loi, ney khong thi hoac tra ve so byte da doc 
}

/*
**Config PCA9685 và set xung PWM để điều khiển servo
*/
static void Init_PCA9685(void)
{
	uint32_t prescaleval, oldmode;
	
	/*
   	Công thức tính prescale value
   	prescale = round (25MHz / (4096*freg)) - 1;
   	trong đó freg là tần số hoạt động truyền vào
   	*/
   	
	I2C_Write(MODE1, 0x00);
	prescaleval = 25000000;
	prescaleval /= 4096;
	prescaleval /= PWM_I2C_Hz;
	prescaleval -= 1.0;
	prescaleval = prescaleval + 3;
	oldmode = I2C_Read(MODE1);
	I2C_Write(MODE1, (oldmode & 0x7F) | 0x10);
	I2C_Write(PRESCALE, prescaleval);
	I2C_Write(MODE1, oldmode);
	msleep(5);
	I2C_Write(MODE1, oldmode | 0x80);
	I2C_Write(MODE2, 0x04);
	PCA9685_setServoPulse(0, 1500);
	PCA9685_setServoPulse(1, 1500);

/* 
Dưới đây là chương trình thay đổi góc quay 1 servo, tăng dần 10 độ mỗi lần di chuyển
yêu cầu: Sinh viên áp dụng để thiết kế chương trình điều khiển cùng lúc 2 servo
Tăng dần 5 độ mỗi lần di chuyển
*/
	while (1)
	{
		int i;
		for (i = 10; i < 90; i = i + 10)
		{
			PCA9685_Set_Rotation_Angle(1, i);
			if (i < 85)
				PCA9685_Set_Rotation_Angle(0, i);
			msleep(500);
		}
	}
	pr_info("Cross Complie");
}


/* 
Hàm set duty cycle( điều chế độ rộng xung )
Tham số truyền vào: - channel: kênh cần gửi tín hiệu trong PCA9685
					   - on: Led ON
					   - off: Led OFF
Lý thuyết: 
	Thời gian active mức 1 của LED và chu kỳ làm việc của PWM có thể được điều
	khiển độc lập bằng cách sử dụng các thanh ghi Led on và led off
	Sẽ có 2 thanh ghi 12 bit cho mỗi đầu ra LED, cả 2 thanh ghi này sẽ lưu giá trị từ 0 - 4095
	một thanh ghi giữ thời gian bật và một thanh ghi giữ thời gian tắt. thời gian bật tắt được so
	sánh với giá trị của bộ đếm 12 bit chạy từ 0000h - 0fffh
*/
static void PCA9685_setPWM(uint8_t channel, uint16_t on, uint16_t off)
{
	I2C_Write(LED0_ON_L + 4 * channel, on & 0xFF);
	I2C_Write(LED0_ON_H + 4 * channel, on >> 8);
	I2C_Write(LED0_OFF_L + 4 * channel, off & 0xFF);
	I2C_Write(LED0_OFF_H + 4 * channel, off >> 8);
}

static void PCA9685_setServoPulse(uint8_t channel, uint16_t value)
{
	value = value * 4096 / 20000;
	PCA9685_setPWM(channel, 0, value);
};

/* 
Hàm set góc quay cho servo
Tham số truyền vào: - channel: Kênh cần gửi tín hiệu
					   - angle: góc quay( từ 10 -> 170 độ)
*/
static void PCA9685_Set_Rotation_Angle(uint8_t channel, uint8_t Angle)
{
	uint16_t temp;
	if (Angle < 0 && Angle > 180)
	{
		pr_info("Angle out of range \n");
	}
	else
	{
		temp = Angle * (2000 / 180) + 500;
		PCA9685_setServoPulse(channel, temp);
	}
}

// Ham nay duoc goi khi Slave duoc tim thay
// Chu y: Chi duoc goi duy nhat mot lan khi load driver.

static int PCA9685_probe(struct i2c_client *client,
                         const struct i2c_device_id *id)
{

	Init_PCA9685();
	pr_info("PCA9685 Init!!!");
	return 0;

}

// Ham nay duoc goi khi salve bi remove
// Chu y: Chi duoc goi duy nhat mot lan khi unload driver.
static int PCA9685_remove(struct i2c_client *client)
{
	PCA9685_Set_Rotation_Angle(0, 80);
	PCA9685_Set_Rotation_Angle(1, 80);
	pr_info("PCA9685 Driver Remove!!!");
	return 0;
}


/*
I2C Board Info structure
*/
static struct i2c_board_info PCA9685_i2c_board_info = {
	I2C_BOARD_INFO(SLAVE_DEVICE_NAME, PCA9685_SLAVE_ADDR)
};

/*
Struct chua salve ID
*/
static const struct i2c_device_id PCA9685_id[] = {
	{ SLAVE_DEVICE_NAME, 0},
	{ }
};
MODULE_DEVICE_TABLE(i2c, PCA9685_id);

/*
I2C driver Structure can add vao linux
*/
static struct i2c_driver PCA9685_driver = {
	.driver = {
		.name = SLAVE_DEVICE_NAME,
		.owner = THIS_MODULE,
	},
	.probe = PCA9685_probe,
	.remove = PCA9685_remove,
	.id_table = PCA9685_id,
};



/*
Module Init function
*/
static int __init PCA9685_driver_init(void)
{
	int ret = -1;
	PCA9685_i2c_adapter = i2c_get_adapter(I2C_BUS_AVAILABLE);

	if (PCA9685_i2c_adapter != NULL)
	{
		PCA9685_i2c_client = i2c_new_client_device(PCA9685_i2c_adapter, &PCA9685_i2c_board_info);

		if (PCA9685_i2c_client != NULL)
		{
			i2c_add_driver(&PCA9685_driver);
			ret = 0;
		}

		i2c_put_adapter(PCA9685_i2c_adapter);
	}

	pr_info("Driver Added!!!\n");
	return ret;
}

/*
Module Exit function
*/
static void __exit PCA9685_driver_exit(void)
{
	i2c_unregister_device(PCA9685_i2c_client);
	i2c_del_driver(&PCA9685_driver);
	pr_info("Driver Removed!!!\n");
}

module_init(PCA9685_driver_init);
module_exit(PCA9685_driver_exit);

MODULE_LICENSE("GPL");
MODULE_AUTHOR("LAB411 <lab411@sis.hust.edu.vn>");
MODULE_DESCRIPTION("I2C driver PCA9685");
MODULE_VERSION("1.34");
