/*******************************************************************************
* ASTex:                                                                       *
* Copyright (C) IGG Group, ICube, University of Strasbourg, France             *
*                                                                              *
* This library is free software; you can redistribute it and/or modify it      *
* under the terms of the GNU Lesser General Public License as published by the *
* Free Software Foundation; either version 2.1 of the License, or (at your     *
* option) any later version.                                                   *
*                                                                              *
* This library is distributed in the hope that it will be useful, but WITHOUT  *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or        *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License  *
* for more details.                                                            *
*                                                                              *
* You should have received a copy of the GNU Lesser General Public License     *
* along with this library; if not, write to the Free Software Foundation,      *
* Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.           *
*                                                                              *
* Web site: https://astex-icube.github.io                                      *
* Contact information: astex@icube.unistra.fr                                  *
*                                                                              *
*******************************************************************************/

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

#include<ASTex/image_common.h>

// forward (must be define whehe .h included
void app_mouse_clicked(int button, int x, int y,int id);
void app_key_pressed(int code, char key, int);

class ImageViewer : public QMainWindow
{
	Q_OBJECT

public:
	/**
	 * @brief Constructor
	 * @param name windows name
	 * @param app appliction ptr (for out with ESC)
	 * @param id id of win (for use of app_mouse_clicked & app_key_pressed)
	 */

	inline ImageViewer(const std::string& name = "") :
		m_imageLabel(new QLabel), title_(name.c_str()), zoom_(1), scale_win_(1), x_(0), y_(0),
		app_mouse_clicked_([](int, int, int) {}),
		app_key_pressed_([](int, char) {})
	{
		m_imageLabel->setBackgroundRole(QPalette::Base);
		m_imageLabel->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
		m_imageLabel->setScaledContents(false);
		setCentralWidget(m_imageLabel);
		setWindowTitle(title_);
	}

	template <typename CB>
	void set_mouse_cb(const CB& cb)
	{
		app_mouse_clicked_ = cb;
	}

	template <typename CB>
	void set_key_cb(const CB& cb)
	{
		app_key_pressed_ = cb;
	}

	template <typename IMG>
	inline auto update(const IMG& img) -> typename std::enable_if<IMG::NB_CHANNELS==1>::type
	{
		using SCALAR = typename IMG::DataType;

		const SCALAR* ptr = img.getDataPtr();
		int w = img.width();
		int h = img.height();

		auto conv = [&] (const SCALAR x) -> uint8_t
		{
			if (std::is_floating_point<SCALAR>::value)
				return ASTex::unnormalized<uint8_t>(x);
			return uint8_t(x);
		};

		QImage im(w,h,QImage::Format_RGB888);
		for(int i=0;i<h;++i)
		{
			uint8_t *optr = im.scanLine(i);
			for(int j=0;j<w;++j)
			{
				*optr++ = conv(*ptr);
				*optr++ = conv(*ptr);
				*optr++ = conv(*ptr++);
			}
		}

		pixmap_ = QPixmap::fromImage(im);
		m_imageLabel->setPixmap(pixmap_);
		m_imageLabel->adjustSize();
		scale_win_ = 1.0;
		while (scale_win_ * pixmap_.height() > 800)
			scale_win_ /= 2;
		scale_window();
	}

	template <typename IMG>
	inline auto update(const IMG& img) -> typename std::enable_if<IMG::NB_CHANNELS==3>::type
	{
		using SCALAR = typename IMG::DataType;

		const SCALAR* ptr = img.getDataPtr();
		int w = img.width();
		int h = img.height();

		auto conv = [&] (const SCALAR x) -> uint8_t
		{
			if (std::is_floating_point<SCALAR>::value)
				return ASTex::unnormalized<uint8_t>(x);
			return uint8_t(x);
		};

		QImage im(w,h,QImage::Format_RGB888);
		for(int i=0;i<h;++i)
		{
			uint8_t *optr = im.scanLine(i);
			for(int j=0;j<w;++j)
			{
				*optr++ = conv(*ptr++);
				*optr++ = conv(*ptr++);
				*optr++ = conv(*ptr++);
			}
		}

		pixmap_ = QPixmap::fromImage(im);
		m_imageLabel->setPixmap(pixmap_);
		m_imageLabel->adjustSize();
		scale_win_ = 1.0;
		while (scale_win_*pixmap_.height() > 800)
			scale_win_ /= 2;
		scale_window();

	}

	template <typename IMG>
	inline auto update(const IMG& img) -> typename std::enable_if<IMG::NB_CHANNELS==4>::type
	{
		using SCALAR = typename IMG::DataType;

		const SCALAR* ptr = img.getDataPtr();
		int w = img.width();
		int h = img.height();

		auto conv = [&] (const SCALAR x) -> uint8_t
		{
			if (std::is_floating_point<SCALAR>::value)
				return ASTex::unnormalized<uint8_t>(x);
			return uint8_t(x);
		};

		QImage im(w,h,QImage::Format_RGB888);
		for(int i=0;i<h;++i)
		{
			uint8_t *optr = im.scanLine(i);
			for(int j=0;j<w;++j)
			{
				*optr++ = conv(*ptr++);
				*optr++ = conv(*ptr++);
				*optr++ = conv(*ptr++);
				*optr++ = conv(*ptr++);
			}
		}
		pixmap_ = QPixmap::fromImage(im);
		m_imageLabel->setPixmap(pixmap_);
		m_imageLabel->adjustSize();
		scale_win_ = 1.0;
		while (scale_win_*pixmap_.height() > 800)
			scale_win_ /= 2;
		scale_window();
	}

	void zoom_update()
	{
		const QPixmap* lpm = m_imageLabel->pixmap();

		int wd = pixmap_.width() / zoom_ * scale_win_;
		int hd = pixmap_.height() / zoom_ * scale_win_;

		real_zoom_x_ = float(pixmap_.width()) / (pixmap_.width() / zoom_);
		real_zoom_y_ = float(pixmap_.width()) / (pixmap_.height() / zoom_);


		x_ = std::max(x_, 0);
		y_ = std::max(y_, 0);
		x_ = std::min(x_, int(pixmap_.width() - wd));
		y_ = std::min(y_, int(pixmap_.height() - hd));
		
		QPixmap npm = pixmap_.copy(x_, y_, wd, hd);
		QPixmap spm = npm.scaled(lpm->width(), lpm->height(), Qt::IgnoreAspectRatio, Qt::FastTransformation);
		m_imageLabel->setPixmap(spm);
		setWindowTitle(title_+ " / zoom="+ QString::number(zoom_)+" x="+ QString::number(x_) + " y="+ QString::number(y_));
		show();
	}

	void scale_window()
	{
		QSize ns = pixmap_.size()*scale_win_;
		setMinimumSize(ns);
		setMaximumSize(ns);
		QPixmap npm(ns);
		m_imageLabel->setPixmap(npm);
		m_imageLabel->adjustSize();
		resize(ns);
		zoom_update();
	}


private:
	QPixmap pixmap_;
	QLabel* m_imageLabel;
	QString title_;
	int zoom_;
	float scale_win_;
	int x_;
	int y_;
	float real_zoom_x_;
	float real_zoom_y_;
	std::function<void(int, int, int)> app_mouse_clicked_;
	std::function<void(int, char)> app_key_pressed_;

	inline void mousePressEvent(QMouseEvent* event)
	{
		app_mouse_clicked_(event->button(),
			x_ + (event->x()) / real_zoom_x_,
			y_ + (event->y() ) / real_zoom_y_);
	}

	inline void keyPressEvent(QKeyEvent* event)
	{
		switch (event->key())
		{
		case Qt::Key_PageUp:
			scale_win_ *= 2;
			zoom_ *= 2;
			scale_window();
			break;
		case Qt::Key_PageDown:
			if (scale_win_ >= 2)
			{
				scale_win_ /= 2;
				zoom_ /= 2;
				scale_window();
			}
			break;
		case Qt::Key_Plus:
			zoom_ *= 2;
			x_ += m_imageLabel->pixmap()->width() / (2 * zoom_);
			y_ += m_imageLabel->pixmap()->height() / (2 * zoom_);
			zoom_update();
			break;
		case Qt::Key_Minus:
			if ((zoom_ >1)&&(zoom_ > scale_win_))
			{
				x_ -= m_imageLabel->pixmap()->width() / (2 * zoom_);
				y_ -= m_imageLabel->pixmap()->height() / (2 * zoom_);
				zoom_ /= 2;
				zoom_update();
			}
			break;
		case Qt::Key_Right:
				++x_;
				zoom_update();
			break;
		case Qt::Key_Left:
				--x_;
				zoom_update();
			break;
		case Qt::Key_Up:
				--y_;
				zoom_update();
				break;
		case Qt::Key_Down:
				++y_;
				zoom_update();
			break;
		case Qt::Key_Escape:
			QApplication::quit();
			break;
		default:
			app_key_pressed_(event->key(), event->text().toStdString()[0]);
		}
	}


	inline void resizeEvent(QResizeEvent* /*event*/)
	{

//		emit(win_resized(event->size()));
	}

//signals:
//   void win_resized(QSize sz);
//
//public slots:
//	inline void win_resize(QSize sz) { resize(sz); }
};




template <typename IMG>
inline std::unique_ptr<ImageViewer> image_viewer(const IMG& img, const std::string& name="", QApplication* app=nullptr)
{
	ImageViewer* view = new ImageViewer(name);
	view->update(img);
	view->show();
	return std::unique_ptr<ImageViewer>(view);
}

#endif	
