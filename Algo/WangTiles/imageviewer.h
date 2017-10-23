
#ifndef IMAGEVIEWER_H
#define IMAGEVIEWER_H

#include <QApplication>
#include <QMainWindow>
#include <QImage>
#include <QtWidgets>
#include <QFileDialog>

#include <cstring>
#include <iostream>
#include <ASTex/internal.h>

QT_BEGIN_NAMESPACE
class QAction;
class QLabel;
class QMenu;
class QScrollArea;
class QScrollBar;
QT_END_NAMESPACE


// forward (must be define whehe .h included
void app_open(const std::string& filename, int id);
void app_mouse_clicked(int button, int x, int y,int id);
void app_key_pressed(int code, char key, int);


class ImageViewer : public QMainWindow
{
	Q_OBJECT

public:
	/**
	 * @brief Constructor
	 * @param name windows name
	 * @param app  ppliction ptr
	 * @param id 
	 */
	ImageViewer(const std::string& name, QApplication* app, int id);
	
	/**
	 * @brief set from Grayu8 image
	 * @param ptr pointer on pixel data
	 * @param w image width
	 * @param h image height
	 * @param zoom int zoom param
	 */
	void set_gray(const unsigned char* ptr,int w,int h, int zoom=1);

	/**
	 * @brief set from RGBu8 image
	 * @param ptr pointer on pixel data
	 * @param w image width
	 * @param h image height
	 * @param zoom int zoom param
	 */
	void set_rgb(const unsigned char* ptr,int w,int h, int zoom=1);

	/**
	 * @brief set from RGBAu8 image
	 * @param ptr pointer on pixel data
	 * @param w image width
	 * @param h image height
	 * @param zoom int zoom param
	 */
	void set_rgba(const unsigned char* ptr,int w,int h, int zoom=1);

	/**
	 * @brief set from Grayf/d image
	 * @param ptr pointer on pixel data
	 * @param w image width
	 * @param h image height
	 * @param zoom int zoom param
	 */
	template <typename REAL>
	void set_gray01(const REAL* ptr,int w,int h, int zoom=1);

	/**
	 * @brief set from Grayf/d image with conversion
	 * @param ptr pointer on pixel data
	 * @param w image width
	 * @param h image height
	 * @param converToRGB lambda (REAL in, uint8_t* out)->void
	 * @param zoom int zoom param
	 */
	template <typename REAL, typename F>
	auto set_gray01(const REAL* ptr,int w,int h, int zoom, const F& converToRGB)
	->  typename std::enable_if<(ASTex::function_traits<F>::arity==2) &&
			(std::is_same<typename ASTex::function_traits<F>::template arg<0>::type,REAL>::value) &&
			(std::is_same<typename ASTex::function_traits<F>::template arg<1>::type,uint8_t*>::value)>::type;

	/**
	* @brief set from RGBd/f image
	* @param ptr pointer on pixel data
	* @param w image width
	* @param h image height
	* @param zoom int zoom param
	*/
	template <typename REAL>
	inline void set_rgb01(const REAL* ptr, int w, int h, int zoom);

	/**
	 * get horizontal scroll bar for easy connect/syncro
	 */
	inline const QScrollBar* hsb() const { return m_scrollArea->horizontalScrollBar(); }
	
	/**
	 * get vertical scroll bar for easy connect/syncro
	 */
	inline const QScrollBar* vsb() const { return m_scrollArea->verticalScrollBar(); }

private:
	QLabel* m_imageLabel;
	QScrollArea* m_scrollArea;
	int m_id;

	inline void mousePressEvent(QMouseEvent* event)
	{
//		app_mouse_clicked(event->button(),
//		                  event->x()-m_scrollArea->x()+m_scrollArea->horizontalScrollBar()->value(),
//		                  event->y()-m_scrollArea->y()+m_scrollArea->verticalScrollBar()->value(),
//		                  m_id);
	}

	inline void keyPressEvent(QKeyEvent* event)
	{
		app_key_pressed(event->key(),event->text().toStdString()[0],m_id);
	}


	inline void resizeEvent(QResizeEvent* event)
	{
		QMainWindow::resizeEvent(event);
		emit(win_resized(event->size()));
	}

signals:
   void win_resized(QSize sz);

public slots:
	inline void win_resize(QSize sz) { resize(sz); }
	void open();
};


inline ImageViewer::ImageViewer(const std::string& name, QApplication* app, int id):
	m_imageLabel(new QLabel), m_scrollArea(new QScrollArea),m_id(id)
{
	m_imageLabel->setBackgroundRole(QPalette::Base);
	m_imageLabel->setSizePolicy(QSizePolicy::Ignored, QSizePolicy::Ignored);
	m_imageLabel->setScaledContents(true);
	m_scrollArea->setBackgroundRole(QPalette::Dark);
	m_scrollArea->setWidget(m_imageLabel);
	m_scrollArea->setVisible(false);
	setCentralWidget(m_scrollArea);
	QMenu* fileMenu = menuBar()->addMenu(tr("&File"));

	QAction *openAct = fileMenu->addAction(tr("&Open"), this, SLOT(open()));
	openAct->setShortcut(QKeySequence::Open);

	QAction *exitAct = fileMenu->addAction(tr("E&xit"), app, SLOT(closeAllWindows()));
	exitAct->setShortcut(tr("Esc"));

	setWindowTitle(QString(name.c_str()));
}


inline void ImageViewer::open()
{
//	QFileDialog dialog(this, tr("Open File"));
//	dialog.selectMimeTypeFilter("image/jpeg");

//	if (dialog.exec() == QDialog::Accepted )
//		app_open(dialog.selectedFiles().first().toStdString(),m_id);
}


inline void ImageViewer::set_gray(const unsigned char* ptr,int w,int h, int zoom)
{
	QImage im(w,h,QImage::Format_RGB888);
	for(int i=0;i<h;++i)
	{
		unsigned char *optr = im.scanLine(i);
		for(int j=0;j<w;++j)
		{
			*optr++ = *ptr;
			*optr++ = *ptr;
			*optr++ = *ptr++;
		}
	}

	if (zoom != 1)
	{
		QImage imz = im.scaledToWidth(w*zoom);
		im.swap(imz);
	}

	m_imageLabel->setPixmap(QPixmap::fromImage(im));
	if ((w<=1024) && (h<=1024))
		m_scrollArea->setMinimumSize(im.width()+2, im.height()+2);
 
	m_scrollArea->setVisible(true);
	m_imageLabel->adjustSize();  
}


inline void ImageViewer::set_rgb(const unsigned char* ptr,int w,int h, int zoom)
{
	QImage im(w,h,QImage::Format_RGB888);
	for(int i=0;i<h;++i)
		std::memcpy(im.scanLine(i), ptr+3*w*i, 3*w);
		
	if (zoom != 1)
	{
		QImage imz = im.scaledToWidth(w*zoom);
		im.swap(imz);
	}

	m_imageLabel->setPixmap(QPixmap::fromImage(im));
	if ((w<=1024) && (h<=1024))
		m_scrollArea->setMinimumSize(im.width()+2, im.height()+2);
 
	m_scrollArea->setVisible(true);
	m_imageLabel->adjustSize();  
}


inline void ImageViewer::set_rgba(const unsigned char* ptr,int w,int h, int zoom)
{
	QImage im(w,h,QImage::Format_RGBA8888);
	for(int i=0;i<h;++i)
		std::memcpy(im.scanLine(i), ptr+3*w*i, 3*w);
		
	if (zoom != 1)
	{
		QImage imz = im.scaledToWidth(w*zoom);
		im.swap(imz);
	}

	m_imageLabel->setPixmap(QPixmap::fromImage(im));
	if ((w<=1024) && (h<=1024))
		m_scrollArea->setMinimumSize(im.width()+2, im.height()+2);
 
	m_scrollArea->setVisible(true);
	m_imageLabel->adjustSize();  
}

template <typename REAL>
inline void ImageViewer::set_gray01(const REAL* ptr,int w,int h, int zoom)
{
	QImage im(w,h,QImage::Format_RGB888);
	for(int i=0;i<h;++i)
	{
		unsigned char *optr = im.scanLine(i);
		for(int j=0;j<w;++j)
		{
			uint8_t val = uint8_t( *ptr * REAL(255));
			*optr++ = val;
			*optr++ = val;
			*optr++ = val;
			ptr++;
		}
	}

	if (zoom != 1)
	{
		QImage imz = im.scaledToWidth(w*zoom);
		im.swap(imz);
	}

	m_imageLabel->setPixmap(QPixmap::fromImage(im));
	if ((w<=1024) && (h<=1024))
		m_scrollArea->setMinimumSize(im.width()+2, im.height()+2);

	m_scrollArea->setVisible(true);
	m_imageLabel->adjustSize();
}


template <typename REAL, typename F>
inline auto ImageViewer:: set_gray01(const REAL* ptr,int w,int h, int zoom, const F& converToRGB)
->  typename std::enable_if<(ASTex::function_traits<F>::arity==2) &&
		(std::is_same<typename ASTex::function_traits<F>::template arg<0>::type,REAL>::value) &&
		(std::is_same<typename ASTex::function_traits<F>::template arg<1>::type,uint8_t*>::value)>::type
{
	QImage im(w,h,QImage::Format_RGB888);
	for(int i=0;i<h;++i)
	{
		unsigned char *optr = im.scanLine(i);
		for(int j=0;j<w;++j)
		{
			converToRGB(*ptr,optr);
			optr+= 3;
			ptr++;
		}
	}

	if (zoom != 1)
	{
		QImage imz = im.scaledToWidth(w*zoom);
		im.swap(imz);
	}
	
	m_imageLabel->setPixmap(QPixmap::fromImage(im));
	if ((w<=1024) && (h<=1024))
		m_scrollArea->setMinimumSize(im.width()+2, im.height()+2);

	m_scrollArea->setVisible(true);
	m_imageLabel->adjustSize();
}


template <typename REAL>
inline void ImageViewer::set_rgb01(const REAL* ptr, int w, int h, int zoom)
{
	QImage im(w, h, QImage::Format_RGB888);
	for (int i = 0; i<h; ++i)
	{
		unsigned char *optr = im.scanLine(i);
		for (int j = 0; j<w; ++j)
		{
			*optr++ = uint8_t(*ptr + *REAL(255));
			*optr++ = uint8_t(*ptr + *REAL(255));
			*optr++ = uint8_t(*ptr + *REAL(255));
		}
	}

	if (zoom != 1)
	{
		QImage imz = im.scaledToWidth(w*zoom);
		im.swap(imz);
	}

	m_imageLabel->setPixmap(QPixmap::fromImage(im));
	if ((w <= 1024) && (h <= 1024))
		m_scrollArea->setMinimumSize(w + 2, h + 2);

	m_scrollArea->setVisible(true);
	m_imageLabel->adjustSize();
}



#endif	
